package org.seqcode.data.readdb;

import java.net.*;
import java.util.*;
import java.io.*;
import java.nio.channels.*;
import javax.security.sasl.*;
import javax.security.auth.callback.*;

/**
 * <p>API for remote access to the readdb server.
 *
 * <p>Calls throw IOException on network errors.
 * Calls throw ClientException on other errors (authentication, authorization, invalid request ,etc.)
 *
 * <p>Client generally assumes that the hit positions are the 5' end of the hit. 
 *
 * <p>Client IS NOT REENTRANT.  Do not overlap calls to a single Client object.
 * 
 * <p>Client connections can be defined as persistent, which keeps the socket open.
 * Use this when expecting many rapid queries (saves connection overhead).
 * By default, connections are opened and closed in the each method call.  
 * 
 * <p>Most method parameters that are object types (eg Integer, Boolean) are optional.  If a null value
 * is passed then no filtering is done based on that parameter.  

 * <p>Standard parameters shared across methods:
 * <ul>
 * <li> alignid is the name of the alignment.
 * <li> isType2 specifies whether to work on ordinary single-ended reads (false) or the "type 2" reads (true). 
 * 		An example of type 2 reads are the R2 reads in paired-end ChIP-exo. 
 * <li> isPaired specifies whether to work on single-ended reads (false) or paired-end reads (true)
 * <li> isLeft if isPaired is true, then isLeft specifies whether to work on the left read (true)
 *     or right read (false) of the pair
 * <li> plusStrand specifies whether to return only reads on the plus strand (true) or minus strand (false). null
 *     means that reads on both strands should be returned.
 * <li> minWeight specifies the minimum weight of reads (or read pairs) to be returned or included in 
 *     the histogram
 * <li> start, stop specify the lowest (inclusive) and highest (exclusive) coordinates of reads
 *     that should be included in the results
 * </ul>
 *
 * 
 */
public class Client implements ReadOnlyClient {    

    public static final String SaslMechanisms[] = {"CRAM-MD5","DIGEST-MD5"};
    Socket socket=null;
    OutputStream outstream;
    BufferedInputStream instream;
    Request request;
    byte[] buffer; //temporary space for receiving data; contents not persistent between method calls
    private static final int BUFFERLEN = 8192*20; //Why is this *20  & server is *16 ?
    private final int socketLoadDataReadTimeout = 1000*60*8; //socket timeout in ms: set to 8 minutes because we should only be relying on the timeout to detect server shutdowns, and some uses of Client (e.g. loading a lot of reads to ReadDB) can take a long time on the Server.
    private final int socketQueryReadTimeout = 60000; //socket timeout in ms for queries
    private boolean printErrors;
    private String hostname, username, passwd;
    private int portnum;
    private boolean persistentConnection=false, chosenPersistentConnection=false;
    private boolean connected=false;
    
    
    /** Connects to a ReadDB server on the specified host and port using the specified 
     * username and password.
     * 
     * Client connections can be defined as persistent, which keeps the socket open.
     * Use this when expecting many rapid queries (saves connection overhead).
     * Don't use this if you expect Client to be left open for a long time without queries 
     * (routers will drop connection).  
     *  
     * @throws IOException on network errors
     * @throws ClientException if the client cannot authenticate to the server
     */
    public Client (String hostname,
                   int portnum,
                   String username,
                   String passwd,
                   boolean persistentConnection) throws IOException, ClientException {
    	this.hostname=hostname;
    	this.portnum = portnum;
    	this.username=username;
    	this.passwd=passwd;
        this.persistentConnection = persistentConnection;
        this.chosenPersistentConnection = persistentConnection;
    }
    public Client (String hostname, int portnum, String username, String passwd) throws IOException, ClientException {
    	this(hostname, portnum, username, passwd, false);
    }
    /**
     * Creates the default connection
     * as specified by ~/.readdb_passwd or a readdb_passwd found in the classpath
     * Must have keys hostname, port, username, and passwd in a format
     * that java.util.PropertyResourceBundle can read
     * 
     * Client connections can be defined as persistent, which keeps the socket open.
     * Use this when expecting many rapid queries (saves connection overhead).
     *
     * @throws IOException on network errors
     * @throws ClientException if the client cannot authenticate to the server
     */
    public Client(boolean persistentConnection) throws IOException, ClientException {
        String homedir = System.getenv("HOME");
        String basename = "readdb_passwd";
        if (System.getenv("READDBROLE") != null) {
            basename = System.getenv("READDBROLE") + basename;
        }
        String fname = homedir + "/." + basename;
        File propfile = new File(fname);
        PropertyResourceBundle bundle = null;
        if (propfile.exists() && propfile.canRead()) {
            bundle = new PropertyResourceBundle(new FileInputStream(propfile));
            if (System.getenv("DEBUGPW") != null) {
                System.err.println("Opening readdb properties from " + propfile);
            }
        } else {
            ClassLoader cl = ClassLoader.getSystemClassLoader();
            URL url = cl.getResource(basename);
            if (System.getenv("DEBUGPW") != null) {
                System.err.println("Opening readdb properties from " + url);
            }            
            if (url != null) {
                bundle = new PropertyResourceBundle(url.openStream());
            } else {
                throw new IOException("Can't read connection properties from " + url);
            }
        }
        this.hostname = bundle.getString("hostname");
        this.portnum = new Integer(bundle.getString("port"));
        this.username = bundle.getString("username");
        this.passwd = bundle.getString("passwd");
        this.persistentConnection = persistentConnection;
        this.chosenPersistentConnection = persistentConnection;
    }
    public Client() throws IOException, ClientException {this(false);}
    
    public void setPersistentConnection(boolean persistent){
    	this.persistentConnection = persistent;
        this.chosenPersistentConnection = persistent;
        if(!persistent && connected) //If the user has chosen to make this connection non-persistent, hang up any open connections
        	closeConnection();
    }
    
    private void openConnection() throws UnknownHostException, IOException, ClientException {
    	
    	if(!connected){
	    	socket=null;
	    	synchronized(this){
	    		socket = new Socket(hostname,portnum);
		    	socket.setTcpNoDelay(true);
		        socket.setSendBufferSize(BUFFERLEN);
		        socket.setReceiveBufferSize(BUFFERLEN);
		        /* linger = true, time = 0 means that the other side gets a reset if the
		           socket is closed on our end (eg, java exits).  We turn this off just before
		           sending a "bye" to allow for a graceful shutdown.  But the RST in
		           other cases lets the server figure out that we've disappeared
		        */
		        socket.setSoLinger(true,0);
		        socket.setSoTimeout(socketQueryReadTimeout);
		        outstream = socket.getOutputStream();
		        outstream.flush();
		        instream = new BufferedInputStream(socket.getInputStream());
		        buffer = new byte[BUFFERLEN];
		        connected = true;
		        
		        if (!authenticate(hostname,username,passwd)) {
		        	throw new ClientException("Authentication Exception Failed");
		        }
		        request = new Request();
	
		        printErrors = false;
	    	}
    	}
    }
    
    /** 
     * Pings the ReadDB server, 
     * @return true if the server pongs 
     */
    public boolean connectionAlive() throws ClientException{
    	String response="";
    	try{
    		openConnection();
        	synchronized (this){
        		request.type="ping";
    	        sendString(request.toString());
    	        response = readLine();
    	        closeConnection();
        	}
	        if (response.equals("pong")) {
	            return true;
	        } else {
	            return false;
	        }
        }catch(IOException e){
        	//SocketException could be generated by a timeout
        	return false;
        }
    }
    
    /**
     * Return some basic information about the server
     * @return string
     */
    public String getServerInfo(){
    	return("ReadDB\t"+hostname+"\t"+portnum+"\t"+username);
    }
    /**
     * Determines whether the client will print error messages to STDERR.  Useful for debugging 
     * but may produce unwanted screen output.
     */
    public void printErrors(boolean b) {printErrors = b;}
    /**
     * performs the SASL authentication exchange with the server.  currently called by the constructor
     */
    private boolean authenticate(String hostname, String username, String password) throws IOException {
        SaslClient sasl = null;
        try {
            sendString(username + "\n");
            Map<String,String> props = new HashMap<String,String>();
            props.put("Sasl.POLICY_NOPLAINTEXT","true");
            props.put("Sasl.POLICY_NOANONYMOUS","true");
            sasl = Sasl.createSaslClient(SaslMechanisms,
                                         username,
                                         "readdb",
                                         hostname,
                                         props,
                                         new ClientCallbackHandler(username,password));
            if (sasl == null) { return false; }
            byte[] response = (sasl.hasInitialResponse() ? sasl.evaluateChallenge(new byte[0]) : new byte[0]);
            byte continueByte = 1;
            while (continueByte != 0) {
                outstream.write((response.length + "\n").getBytes());
                outstream.write(response);
                outstream.flush();
                int length = Integer.parseInt(readLine());
                byte[] challenge = new byte[length];
                int read = 0;
                while (read < length) {
                    read += instream.read(challenge, read, length - read);
                }
                /* the continueByte tells us whether the server
                   expects to do another round.  Necessary because sometimes
                   isComplete() returned true here but the server wasn't done
                */
                continueByte = (byte)instream.read();
                if (!sasl.isComplete()) { 
                    response = sasl.evaluateChallenge(challenge);
                } else {
                    response = new byte[0];
                }
            }
            sasl.dispose();
            String status = readLine();
            return (status.equals("authenticated as " + username));
        } catch (SaslException e) {
            e.printStackTrace();
            if (sasl != null) {
                sasl.dispose();
            }
            return false;
        } catch (ClientException e) {
			e.printStackTrace();
			return false;
		}
    }
    
    /** sends a string to the server and flushes the socket 
     */
    private void sendString(String s) throws IOException, ClientException {
    	if(!connected)
    		throw new ClientException("Client not connected!");
    	outstream.write(s.getBytes());
        outstream.flush();
    }
    /** reads one line from the server.  blocking.
     */
    private String readLine() throws IOException {
        int pos = 0;
        int i;
        while ((i = instream.read()) != -1) {
            if (i == '\n') {
                break;
            } else {
                buffer[pos++] = (byte)i;
            }
        }
        String out = new String(buffer,0,pos);
        return out;
    }
    /**
     * Tells the server to shut itself down.  Use this to stop the server process.
     * @throws IOException on network errors
     * @throws ClientException if the user isn't authorized to shut the server down.
     */
    public void shutdown() throws IOException, ClientException{
    	openConnection();
    	synchronized(this){
    		request.type = "shutdown";
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            System.out.println(response);
	        }else{
	        	System.out.println("ReadDB shutting down");
	        }
    	}
    	closeConnection();
    }
    /**
     * Tells the server to clear its cache. Only used for debugging & testing.
     * @throws IOException on network errors
     * @throws ClientException if the user isn't authorized to shut the server down.
     */
    public void clearServerCache() throws IOException, ClientException{
    	openConnection();
    	synchronized(this){
    		request.type = "clearcache";
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            System.out.println(response);
	        }else{
	        	System.out.println("ReadDB cache cleared");
	        }
    	}
    	closeConnection();
    }
    /** this was to fix a bug in the server.  You shouldn't need it for general use.
     * Regenerate the index for this alignment and chromosome
     */
    public void reIndex(String align, int chrom, boolean isType2, boolean paired) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type = "reindex";
	        request.alignid = align;
	        request.chromid = chrom;
	        request.isType2 = isType2;
	        request.isPaired = paired;
	        sendString(request.toString());
	        outstream.flush();
	        String response = readLine();
	        if (!response.equals("OK")) {
	            System.out.println(response);
	        }
    	}
    	closeConnection();
    }
    /** This was to fix a bug in the server.  You shouldn't need it for general use.
     * Resort the hits for a single-ended alignment and regenerate the index.
     */
    public void checksort(String align, int chrom) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.clear();
	        request.type = "checksort";
	        request.alignid = align;
	        request.chromid = chrom;
	        sendString(request.toString());
	        outstream.flush();
	        String response = readLine();
	        if (!response.equals("OK")) {
	            System.out.println(response);
	        }
    	}
    	closeConnection();
    }
    /**
     * Stores a set of SingleHit objects (representing an un-paired or single-ended read
     * aligned to a genome) in the specified alignment.  The hits are appended
     * to any hits that have already been stored in the alignment.
     
     */
    public void storeSingle(String alignid, List<SingleHit> allhits, boolean isType2) throws IOException, ClientException {
    	openConnection();
    	socket.setSoTimeout(socketLoadDataReadTimeout);

    	WritableByteChannel outchannel = Channels.newChannel(outstream);
	    int step = 10000000;
        for (int pos = 0; pos < allhits.size(); pos += step) {
            Map<Integer, List<SingleHit>> map = new HashMap<Integer,List<SingleHit>>();
            for (int i = pos; i < pos + step && i < allhits.size(); i++) {
                SingleHit h = allhits.get(i);
                if (!map.containsKey(h.chrom)) {
                    map.put(h.chrom, new ArrayList<SingleHit>());
                }
                map.get(h.chrom).add(h);
            }
            for (int chromid : map.keySet()) {
            	synchronized(this){
            	    List<SingleHit> hits = map.get(chromid);
	                Collections.sort(hits);
	                int chunk = step;
	                for (int startindex = 0; startindex < hits.size(); startindex += chunk) {
	                    int count = ((startindex + chunk) < hits.size()) ? chunk : (hits.size() - startindex);
	                    request.clear();
	                    request.type="storesingle";
	                    request.alignid=alignid;
	                    request.chromid = chromid;
	                    request.isType2 = isType2;
	                    request.map.put("numhits",Integer.toString(count));
	                    try{
	                    	sendString(request.toString());
		                    String response = readLine();
		                    if (!response.equals("OK")) {
		                        System.err.println("not-OK response to request: " + response);
		                        System.err.println("request was " + request);
		                        throw new ClientException(response);
		                    }
		                    int[] ints = new int[count];
		                    for (int i = startindex; i < startindex + count; i++) {
		                        ints[i - startindex] = hits.get(i).pos;
		                    }   
		                    Bits.sendIntsByChannel(ints, outchannel);
		                    float[] floats = new float[count];
		                    for (int i = startindex; i < startindex + count; i++) {
		                        floats[i - startindex] = hits.get(i).weight;
		                        ints[i - startindex] = Hits.makeLAS(hits.get(i).length, hits.get(i).strand);
		                    }
		                    Bits.sendFloatsByChannel(floats, outchannel);
		                    Bits.sendIntsByChannel(ints, outchannel);
		            
		                    System.err.println("Sent " + count + " hits to the server for " + chromid + "," + alignid);
		                    outstream.flush();        
		                    response = readLine();
		                    if (!response.equals("OK")) {
		                        throw new ClientException(response);
		                    }
	                    }catch (IOException ioe){
	                    	//IOException here is probably a socket time-out. 
	                    	//I think it's best to kill the process at this point, since we won't know if the sent reads actually got loaded.
	                    	System.err.println(ioe);
	                    	System.exit(1);
	                    }
	                }
	            }
	        }
    	}
        socket.setSoTimeout(socketQueryReadTimeout);
        
        closeConnection();
    }
    /**
     * Stores a set of PairedHit objects (representing an paired-ended read
     * aligned to a genome) in the specified alignment.  The hits are appended
     * to any hits that have already been stored in the alignment
     */
    public void storePaired(String alignid, List<PairedHit> allhits) throws IOException, ClientException {
    	openConnection();
    	socket.setSoTimeout(socketLoadDataReadTimeout);
    	
    	WritableByteChannel outchannel = Channels.newChannel(outstream);
    	Map<Integer, List<PairedHit>> map = new HashMap<Integer,List<PairedHit>>();
    	for (PairedHit h : allhits) {
            if (!map.containsKey(h.leftChrom)) {
                map.put(h.leftChrom, new ArrayList<PairedHit>());
            }
            map.get(h.leftChrom).add(h);
        }
        for (int chromid : map.keySet()) {
        	synchronized(this){
	            List<PairedHit> hits = map.get(chromid);
	            System.err.println("SENDING PAIRED HITS n="+hits.size() + " for chrom " + chromid);
	            int chunk = 1000000;
	            for (int startindex = 0; startindex < hits.size(); startindex += chunk) {
	                int count = ((startindex + chunk) < hits.size()) ? chunk : (hits.size() - startindex);
	
	                request.clear();
	                request.type="storepaired";
	                request.alignid=alignid;
	                request.chromid=chromid;
	                request.isLeft=true;
	                request.map.put("numhits",Integer.toString(count));
	                try{
		                sendString(request.toString());
		                String response = readLine();
		                if (!response.equals("OK")) {
		                    System.err.println("not-OK response to request: " + response);
		                    System.err.println("request was " + request);
		                    throw new ClientException(response);
		                }
		                int[] ints = new int[count];
		                for (int i = startindex; i < startindex + count; i++) {
		                    ints[i-startindex] = hits.get(i).leftPos;
		                }        
		                Bits.sendIntsByChannel(ints, outchannel);
		                float[] floats = new float[count];
		                int[] codes = new int[count];
		                for (int i = startindex; i < startindex + count; i++) {
		                	floats[i-startindex] = hits.get(i).weight;
		                    codes[i-startindex] = hits.get(i).pairCode;
		                    ints[i-startindex] = Hits.makeLAS(hits.get(i).leftLength, hits.get(i).leftStrand,
		                                                      hits.get(i).rightLength, hits.get(i).rightStrand);
		
		                }
		                Bits.sendFloatsByChannel(floats, outchannel);
		                Bits.sendIntsByChannel(codes, outchannel);
		                Bits.sendIntsByChannel(ints, outchannel);
		                for (int i = startindex; i < startindex + count; i++) {
		                    ints[i-startindex] = hits.get(i).rightChrom;
		                }        
		                Bits.sendIntsByChannel(ints, outchannel);
		                for (int i = startindex; i < startindex + count; i++) {
		                    ints[i-startindex] = hits.get(i).rightPos;
		                }        
		                Bits.sendIntsByChannel(ints, outchannel);
		                System.err.println("Sent " + count + " hits to the server");
		                outstream.flush();        
		                response = readLine();
		                if (!response.equals("OK")) {
		                    if (printErrors) {
		                        System.err.println("not-OK response to request: " + response);
		                        System.err.println("request was " + request);
		                    }
		                    throw new ClientException(response);
		                }
	                }catch (IOException ioe){
                    	//IOException here is probably a socket time-out. 
                    	//I think it's best to kill the process at this point, since we won't know if the sent reads actually got loaded.
                    	System.err.println(ioe);
                    	System.exit(1);
                    }
	            }
	        }
    	}
        socket.setSoTimeout(socketQueryReadTimeout);
        
        closeConnection();
    }
    
    /** Returns true if the alignment and chromosome exist and are accessible
     * to the user.  Returns false if they don't exist or if they aren't accessible
     */
    public boolean exists(String alignid) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.clear(); 
	        request.type="exists";
	        request.alignid=alignid;
	        sendString(request.toString());
	        String response = readLine();
	        closeConnection();
	        if (response.equals("exists")) {
	            return true;
	        } else if (response.equals("unknown")) {
	            return false;
	        } else {
	            return false;
	        }
    	}
    }
    
    /**
     * Deletes an alignment (all chromosomes).  isPaired specifies whether to delete
     * the paired or single ended reads.
     */
    public void deleteAlignment(String alignid, boolean isPaired) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="deletealign";
	        request.isPaired = isPaired;
	        request.alignid=alignid;
	        sendString(request.toString());
	        String response = readLine();
	        closeConnection();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
    	}
    }
    /**
     * Returns the set of chromosomes that exist for this alignment. 
     */
    public Set<Integer> getChroms(String alignid, boolean isType2, boolean isPaired, Boolean isLeft) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="getchroms";
	        request.isType2 = isType2;
	        request.isLeft = isLeft;
	        request.isPaired = isPaired;
	        request.alignid=alignid;
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        int numchroms = Integer.parseInt(readLine());
	        Set<Integer> output = new HashSet<Integer>();
	        while (numchroms-- > 0)
	            output.add(Integer.parseInt(readLine()));
	        closeConnection();
	        return output;
    	}
    }
    /**
     * Returns the total number of hits in this alignment.  
     */
    public int getCount(String alignid, boolean isType2, boolean isPaired, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
    	persistentConnection=true; //Override to keep Socket open for faster multi-queries
    	int count = 0;
        for (int c : getChroms(alignid, isType2, isPaired, isLeft)) {
            count += getCount(alignid, c, isType2, isPaired, null,null,null,isLeft,plusStrand);
        }
        persistentConnection=chosenPersistentConnection;
        closeConnection();
        return count;
    }
    /**
     * Returns the sum of the weights of all hits in this alignment
     */
    public double getWeight(String alignid,  boolean isType2, boolean isPaired, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
    	persistentConnection=true; //Override to keep Socket open for faster multi-queries
    	double total = 0;
        for (int c : getChroms(alignid, isType2, isPaired, isLeft)) {
            total += getWeight(alignid, c, isType2, isPaired, null, null, null, isLeft, plusStrand);
        }
        persistentConnection=chosenPersistentConnection;        
        closeConnection();
        return total;
    }
    /**
     * Returns the total number of unique positions in this alignment.  
     */
    public int getNumPositions(String alignid, boolean isType2, boolean isPaired, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
    	persistentConnection=true; //Override to keep Socket open for faster multi-queries
    	int pos = 0;
        for (int c : getChroms(alignid, isType2, isPaired, isLeft)) {
            pos += getNumPositions(alignid, c, isType2, isPaired, null,null,null,isLeft,plusStrand);
        }
        persistentConnection=chosenPersistentConnection;
        closeConnection();
        return pos;
    }
    /**
     * Returns the total number of unique paired positions in this alignment.  
     */
    public int getNumPairedPositions(String alignid, boolean isType2, Boolean isLeft) throws IOException, ClientException {
    	persistentConnection=true; //Override to keep Socket open for faster multi-queries
    	int pos = 0;
        for (int c : getChroms(alignid, isType2, true, isLeft)) {
            pos += getNumPairedPositions(alignid, c, isType2, null,null,null,isLeft);
        }
        persistentConnection=chosenPersistentConnection;
        closeConnection();
        return pos;
    }
    

    /** returns the total number of hits on the specified chromosome in the alignment.
     * Any of the object parameters can be set to null to specify "no value"
     */
    public int getCount(String alignid, int chromid, boolean isType2, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand)  throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="count";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isPaired = paired;
	        request.isLeft = isLeft == null ? true : isLeft;
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        int numhits = Integer.parseInt(readLine());
	        closeConnection();
	        return numhits;
    	}
    }
    /** returns the total weight on the specified chromosome in this alignment
     */
    public double getWeight(String alignid, int chromid, boolean isType2, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="weight";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isPaired = paired;
	        request.isLeft = isLeft == null ? true : isLeft;
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        double weight = Double.parseDouble(readLine());
	        closeConnection();
	        return weight;
    	}
    }
    /** returns the total number of unique positions on the specified chromosome in the alignment.
     * Any of the object parameters can be set to null to specify "no value"
     */
    public int getNumPositions(String alignid, int chromid, boolean isType2, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand)  throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="numpositions";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isPaired = paired;
	        request.isLeft = isLeft == null ? true : isLeft;
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        int numpos = Integer.parseInt(readLine());
	        closeConnection();
	        return numpos;
    	}
    }
    /** returns the total number of unique paired positions on the specified chromosome in the alignment.
     * Note that total number of unique paired positions will depend on the setting of "isLeft": 
     *    e.g. will count pair positions with left read in specified chromosome.
     * Any of the object parameters can be set to null to specify "no value"
     */
    public int getNumPairedPositions(String alignid, int chromid, boolean isType2, Integer start, Integer stop, Float minWeight, Boolean isLeft)  throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="numpairpositions";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = null;
	        request.isPaired = true;
	        request.isLeft = isLeft == null ? true : isLeft;
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        int numpos = Integer.parseInt(readLine());
	        closeConnection();
	        return numpos;
    	}
    }
    /** 
     * returns the sorted (ascending order) hit positions in the specified range of a chromosome,alignment pair.
     */ 
    public int[] getPositions(String alignid, int chromid, boolean isType2, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="gethits";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isPaired = paired;
	        request.isLeft = isLeft;
	        request.map.put("wantpositions","1");
	        sendString(request.toString());        
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        int numhits = Integer.parseInt(readLine());
	        IntBP ints = new IntBP(numhits);
	        ReadableByteChannel rbc = Channels.newChannel(instream);
	        Bits.readBytes(ints.bb, rbc);
	        int[] pos = new int[numhits];
	        for (int i = 0; i < numhits; i++) {
	        	pos[i]=ints.get(i);
	        }
	    	
	        closeConnection();
	        return pos;
    	}
    }
    /** 
     * returns the hit weights in the specified range of a chromosome,alignment pair.  The weights
     * will be in the same order as the sorted positions returned by getPositions()
     */ 
    public float[] getWeightsRange(String alignid, int chromid, boolean isType2, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="gethits";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isPaired = paired;
	        request.isLeft = isLeft;
	        request.map.put("wantweights","1");
	        sendString(request.toString());        
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        int numhits = Integer.parseInt(readLine());
	        FloatBP flts = new FloatBP(numhits);
	        ReadableByteChannel rbc = Channels.newChannel(instream);
	        Bits.readBytes(flts.bb, rbc);
	        float[] wr = new float[numhits];
	        for (int i = 0; i < numhits; i++) {
	        	wr[i]=flts.get(i);
	        }
	        
	        closeConnection();
	        return wr;
    	}
    }
    public List<SingleHit> getSingleHits(String alignid, int chromid, boolean isType2, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="gethits";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isPaired = false;
	        request.map.put("wantpositions","1");
	        request.map.put("wantweights","1");
	        request.map.put("wantlengthsandstrands","1");
	        sendString(request.toString());        
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        List<SingleHit> output = new ArrayList<SingleHit>();
	        int numhits = Integer.parseInt(readLine());
	        for (int i = 0; i < numhits; i++) {
	            output.add(new SingleHit(chromid,0,(float)0.0,false,(short)0));
	        }
	        IntBP ints = new IntBP(numhits);
	        ReadableByteChannel rbc = Channels.newChannel(instream);
	        ints.bb.clear();
	        Bits.readBytes(ints.bb, rbc);
	        for (int i = 0; i < numhits; i++) {
	            output.get(i).pos = ints.get(i);
	        }
	        FloatBP floats = new FloatBP(numhits);
	        floats.bb.clear();
	        Bits.readBytes(floats.bb, rbc);
	        for (int i = 0; i < numhits; i++) {
	            output.get(i).weight = floats.get(i);
	        }
	        ints.bb.clear();
	        Bits.readBytes(ints.bb, rbc);
	        for (int i = 0; i < numhits; i++) {
	            int j = ints.get(i);
	            SingleHit h = output.get(i);
	            h.length = Hits.getLengthOne(j);
	            h.strand = Hits.getStrandOne(j);
	        }
	        closeConnection();
	        return output;
    	}
    }
    public List<PairedHit> getPairedHits(String alignid, int chromid, boolean isLeft, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="gethits";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isPlusStrand = plusStrand;
	        request.isLeft = isLeft;
	        request.isType2 = false;
	        request.isPaired = true;
	        request.map.put("wantpositions","1");
	        request.map.put("wantweights","1");
	        request.map.put("wantpaircodes","1");
	        request.map.put("wantlengthsandstrands","1");
	        request.map.put("wantotherchroms","1");
	        request.map.put("wantotherpositions","1");
	        sendString(request.toString());        
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(String.format("align %s chrom %d: %s", alignid, chromid, response));
	        }
	        List<PairedHit> output = new ArrayList<PairedHit>();
	        int numhits = Integer.parseInt(readLine());
	        for (int i = 0; i < numhits; i++) {
	            output.add(new PairedHit(chromid,0,false,(short)0,
	                                     chromid,0,false,(short)0,(float)0,0));
	        }
	        IntBP ints = new IntBP(numhits);
	        ReadableByteChannel rbc = Channels.newChannel(instream);
	        Bits.readBytes(ints.bb, rbc);
	        if (isLeft) {
	            for (int i = 0; i < numhits; i++) {
	                output.get(i).leftPos = ints.get(i);
	            }
	        } else {
	            for (int i = 0; i < numhits; i++) {
	                output.get(i).rightPos = ints.get(i);
	            }
	        }
	        FloatBP floats = new FloatBP(numhits);
	        Bits.readBytes(floats.bb, rbc);
	        for (int i = 0; i < numhits; i++) {
	            output.get(i).weight = floats.get(i);
	        }
	
	        Bits.readBytes(ints.bb, rbc);
	        for (int i = 0; i < numhits; i++) {
	            output.get(i).pairCode = ints.get(i);
	        }
	        
	        Bits.readBytes(ints.bb, rbc);
	        if (isLeft) {
	            for (int i = 0; i < numhits; i++) {
	                int j = ints.get(i);
	                PairedHit h = output.get(i);                
	                h.leftLength = Hits.getLengthOne(j);
	                h.leftStrand = Hits.getStrandOne(j);
	                h.rightLength = Hits.getLengthTwo(j);
	                h.rightStrand = Hits.getStrandTwo(j);
	            }
	        } else {
	            for (int i = 0; i < numhits; i++) {
	                int j = ints.get(i);
	                PairedHit h = output.get(i);                
	                h.leftLength = Hits.getLengthTwo(j);
	                h.leftStrand = Hits.getStrandTwo(j);
	                h.rightLength = Hits.getLengthOne(j);
	                h.rightStrand = Hits.getStrandOne(j);
	            }
	        }
	        Bits.readBytes(ints.bb, rbc);
	        if (isLeft) {
	            for (int i = 0; i < numhits; i++) {
	                output.get(i).rightChrom = ints.get(i);
	            }
	        } else {
	            for (int i = 0; i < numhits; i++) {
	                output.get(i).leftChrom = ints.get(i);
	            }
	        }
	
	        Bits.readBytes(ints.bb, rbc);
	        if (isLeft) {
	            for (int i = 0; i < numhits; i++) {
	                output.get(i).rightPos = ints.get(i);
	            }
	        } else {
	            for (int i = 0; i < numhits; i++) {
	                output.get(i).leftPos = ints.get(i);
	            }
	        }
	        closeConnection();
	        return output;
    	}
    }


    /**
     * returns a TreeMap from positions to counts representing a histogram
     * of the hits in a range with the specified binsize.  Bins with a count
     * of zero are not included in the output.
     * 
     * Ex: getHistgram("foo","1+",1,100,10)
     *     6: 5
     *    16: 0
     *    36: 30
     *
     * minweight is the minimum weight for reads to be included in the histogram.
     *
     * dedup is the limit on how many times reads with any given 5' position will be counted.
     *  A value of zero means no limit.  A limit of, eg, 2, means that at most two reads at
     *  any 5' position will be included in the output.  For weighted histograms, the choice of reads
     *  included is unspecified.  For methods that operate on a set of alignments, this many
     *  reads from each alignment will be included.
     *
     * Normally, a read is only counted in a single bin as defined by its position (generally
     * the 5' end of the read).  A non-zero read-Extension counts the read in any
     * bin that you cross within that many bases of the read's position.  A negative value
     * goes backwards (smaller coordinates) and a larger value goes forward.  You need
     * to get the sign right depending on the strandedness of the chromosome that you're
     * working with.
     */
    public TreeMap<Integer,Integer> getHistogram(String alignid, int chromid, boolean isType2, boolean paired, int extension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getHistogram(alignid, chromid, isType2, paired, extension,binsize,0,start,stop,minWeight,plusStrand,true);
    }
    public TreeMap<Integer,Integer> getHistogram(String alignid, int chromid, boolean isType2, boolean paired, int extension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand, boolean isLeft) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="histogram";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.isLeft = isLeft;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isPaired = paired;
	        request.map.put("binsize",Integer.toString(binsize));
	        if (dedup > 0) {
	            request.map.put("dedup",Integer.toString(dedup));
	        }
	        if (extension != 0) {
	            request.map.put("extension",Integer.toString(extension));
	        }
	        sendString(request.toString());        
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        TreeMap<Integer,Integer> output = new TreeMap<Integer,Integer>();
	        int numints = Integer.parseInt(readLine());
	        if(numints>0){
		        IntBP ints = new IntBP(numints);
		        ReadableByteChannel rbc = Channels.newChannel(instream);
		        Bits.readBytes(ints.bb, rbc);
		        for (int i = 0; i < numints; i+=2) {
		        	output.put(ints.get(i), ints.get(i+1));
		        }
	        }
	        closeConnection();
	        return output;
    	}
    }
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, int chromid, boolean isType2, boolean paired, int extension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getWeightHistogram(alignid, chromid, isType2, paired, extension, binsize, 0, start,stop,minWeight,plusStrand, true);
    }
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, int chromid, boolean isType2, boolean paired, int extension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand, boolean isLeft) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="weighthistogram";
	        request.alignid=alignid;
	        request.chromid=chromid;
	        request.start = start;
	        request.end = stop;
	        request.minWeight = minWeight;
	        request.isType2 = isType2;
	        request.isPlusStrand = plusStrand;
	        request.isLeft = isLeft;
	        request.isPaired = paired;
	        request.map.put("binsize",Integer.toString(binsize));
	        if (dedup > 0)
	            request.map.put("dedup",Integer.toString(dedup));
	        if (extension!=0)
	            request.map.put("extension",Integer.toString(extension));
	        
	        sendString(request.toString());        
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        TreeMap<Integer,Float> output = new TreeMap<Integer,Float>();
	        int numints = Integer.parseInt(readLine());
	        if(numints>0){
	        	ReadableByteChannel rbc = Channels.newChannel(instream);
	        	IntBP ints = new IntBP(numints);
		        Bits.readBytes(ints.bb, rbc);
		        FloatBP floats = new FloatBP(numints);
		        Bits.readBytes(floats.bb, rbc);
		        for (int i = 0; i < numints; i++) {
		        	output.put(ints.get(i), floats.get(i));
		        }
	        }
	        closeConnection();
	        return output;
    	}
    }

    public TreeMap<Integer,Integer> getHistogram(Collection<String> alignids, int chromid, boolean isType2, boolean paired, int extension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getHistogram(alignids,chromid,isType2, paired,extension,binsize,0,start,stop,minWeight,plusStrand);
    }
    public TreeMap<Integer,Integer> getHistogram(Collection<String> alignids, int chromid, boolean isType2, boolean paired, int extension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
    	persistentConnection=true; //Override to keep Socket open for faster multi-queries
    	TreeMap<Integer,Integer> output = null;
        for (String alignid : alignids) {
            TreeMap<Integer,Integer> o = getHistogram(alignid,chromid,isType2, paired,extension,binsize,dedup,start,stop,minWeight,plusStrand,true);
            if(paired) //run for isLeft =true & false
            	o.putAll(getHistogram(alignid,chromid,isType2, paired,extension,binsize,dedup,start,stop,minWeight,plusStrand,false));
            for (int k : o.keySet()) { 
                if ((k - start - binsize / 2) % binsize != 0 ) {
                    System.err.println(String.format("Bad key %d for binsize %d and start %d in %s,%d",
                                                     k, binsize, start, alignid,chromid));
                }
            }

            if (output == null) {
                output = o;
            } else {
                for (int k : o.keySet()) {
                    if (output.containsKey(k)) {
                        output.put(k, output.get(k) + o.get(k));
                    } else {
                        output.put(k,o.get(k));
                    }
                }
            }    
        }
        persistentConnection=chosenPersistentConnection;
        closeConnection();
        return output;
    }
    public TreeMap<Integer,Float> getWeightHistogram(Collection<String> alignids, int chromid, boolean isType2, boolean paired, int extension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getWeightHistogram(alignids,chromid,isType2, paired,extension,binsize,0,start,stop,minWeight,plusStrand);
    }
    public TreeMap<Integer,Float> getWeightHistogram(Collection<String> alignids, int chromid, boolean isType2, boolean paired, int extension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
    	persistentConnection=true; //Override to keep Socket open for faster multi-queries
    	TreeMap<Integer,Float> output = null;
        for (String alignid : alignids) {
            TreeMap<Integer,Float> o = getWeightHistogram(alignid,chromid,isType2, paired,extension,binsize,dedup,start,stop,minWeight,plusStrand, true);
            if(paired) //run for isLeft =true & false
            	o.putAll(getWeightHistogram(alignid,chromid,isType2, paired,extension,binsize,dedup,start,stop,minWeight,plusStrand, false));
            if (output == null) {
                output = o;
            } else {
                for (int k : o.keySet()) {
                    if (output.containsKey(k)) {
                        output.put(k, output.get(k) + o.get(k));
                    } else {
                        output.put(k,o.get(k));
                    }
                }
            }            
        }
        persistentConnection=chosenPersistentConnection;
        closeConnection();
        return output;
    }

    /**
     * Returns a Map from READ, WRITE, and ADMIN to lists of principals that have those privileges on the specified alignment.
     */
    public Map<String,Set<String>> getACL(String alignid) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="getacl";
	        request.alignid=alignid;
	        sendString(request.toString());
	        String response = readLine();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
	        Map<String,Set<String>> output = new HashMap<String,Set<String>>();
	        fillPartACL(output);
	        fillPartACL(output);
	        fillPartACL(output);
	        closeConnection();
	        return output;
    	}
    }
    /* fills one section of the acl output data structure.  A section
       is either read, write, or admin. 
    */
    private void fillPartACL(Map<String,Set<String>> output) throws IOException, ClientException {
    	synchronized(this){
    		String type = readLine();
    		int entries = Integer.parseInt(readLine());
	        Set<String> out = new HashSet<String>();
	        while (entries-- > 0) {
	            out.add(readLine());
	        }
	        output.put(type,out);
    	}
    }
    /**
     * Applies the specified ACLChangeEntry objects to the acl for this experiment/chromosome.
     */
    public void setACL(String alignid, Set<ACLChangeEntry> changes) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="setacl";
	        request.alignid=alignid;
	        for (ACLChangeEntry a : changes) {
	            request.list.add(a.toString());
	        }    
	        sendString(request.toString());
	        String response = readLine();
	        closeConnection();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
    	}
    }
    /**
     * Adds the specified user (princ) to a group.
     */
    public void addToGroup(String princ, String group) throws IOException, ClientException {
    	openConnection();
    	synchronized(this){
    		request.type="addtogroup";
	        request.map.put("princ",princ);
	        request.map.put("group",group);
	        sendString(request.toString());
	        String response = readLine();
	        closeConnection();
	        if (!response.equals("OK")) {
	            if (printErrors) {
	                System.err.println("not-OK response to request: " + response);
	                System.err.println("request was " + request);
	            }
	            throw new ClientException(response);
	        }
    	}
    }
    /**
     * Closes this connection to the server.
     */
    public void close() {
    	if(connected)
        	closeConnection(true);        
    } 
    /**
     * Closes this connection to the server.
     * finalClose is used to force closure in cases where connection is persistant. 
     */
    public void closeConnection(boolean finalClose) {
    	if(!persistentConnection || finalClose){
	    	if (socket == null || !connected) {
	            return;
	        }
	        try {
	        	//Setting linger to false enables an orderly Socket close()
	            socket.setSoLinger(false,0);
	            request.clear();
	            request.type="bye";
	            sendString(request.toString());
	            
	        } catch (IOException e) {
	            e.printStackTrace();
	        } catch (ClientException e) {
				e.printStackTrace();
			}finally{
            	try{
            		outstream.close();
    	            outstream = null;
            	} catch (IOException e) {
    	            e.printStackTrace();
    	        }finally{
            		try {
        	            instream.close();
        	            instream = null;
        	        } catch (IOException e) {
        	            e.printStackTrace();
        	        }finally{
	        	        try {
	        	            socket.close();
	        	            socket = null;
	        	        } catch (IOException e) {
	        	            e.printStackTrace();
	        	        }finally{
	        	        	connected=false;
	        	        }
            		}
            	}
            }
    	}
    } 
    public void closeConnection() {this.closeConnection(false);}

}

/**
 * SASL callback handler for the authentication exchange
 */
class ClientCallbackHandler implements CallbackHandler {
    private String name, pass;
    public ClientCallbackHandler(String n, String p) {name = n; pass = p;}
    public void handle(Callback[] callbacks) {
        for (int i = 0; i < callbacks.length; i++) {
            if (callbacks[i] instanceof NameCallback) {
                NameCallback nc = (NameCallback)callbacks[i];
                nc.setName(name);
            }
            if (callbacks[i] instanceof PasswordCallback) {
                PasswordCallback pc = (PasswordCallback)callbacks[i];
                pc.setPassword(pass.toCharArray());
            }
        }
    }
}