package org.seqcode.data.readdb;

import java.io.*;
import java.nio.channels.*;
import java.net.*;
import java.util.*;
import java.util.logging.*;
import javax.security.sasl.*;
import javax.security.auth.callback.*;

/** 
 * ServerTask represents a client connection.  Server creates ServerTasks when it receives
 * a connection and passes them to Dispatch.  Dispatch manages a pool of WorkerThreads and
 * assigns them to ServerTasks as the tasks appear to be available. 
 * 
 * The current ServerTasks shut themselves down if they haven't been run in a given time period. 
 *   
 */

public class ServerTask {
    
    /* limit in milliseconds before a task connection is closed by server */
    private int taskInactivityLimit; 
    private boolean limitInactivity; 
    /* instance variables */
    /* Server that we're working for.  Need a reference to it so we can ask it
       for paths and configuration information and such
    */
    private Server server;
    /* instance variables set for each request */
    /* true if the client sent a close command.  If the client sent a close command or
       the thread has detected an error from the client, then shouldClose is true
       and the server will close the connection.
    */
    private boolean shouldClose;
    /* Socket, streams from the socket */
    private Socket socket;
    private int haventTriedRead;
    private BufferedInputStream instream;
    private OutputStream outstream;
    private WritableByteChannel outchannel;
    /* if authenticate was successful, this holds a username.  Null otherwise */
    private String username;
    /* buffer for readLine */
    private int bufferpos;
    private byte[] buffer;
    private static final int MAXPARAMLINES = 100;
    /* other variables maintained across calls to Run but reset between connections */
    private Request request;
    private List<String> args;
    private SaslServer sasl;
    private Map<String,String> saslprops;
    private String uname; // temporary, used by authenticate
    private long lastActivity=0;

    public ServerTask(Server serv, Socket s, int inactivityLimit) throws IOException {
        buffer = new byte[8192];
        request = new Request();
        args = new ArrayList<String>();
        saslprops = new HashMap<String,String>();
        saslprops.put("Sasl.POLICY_NOPLAINTEXT","true");
        saslprops.put("Sasl.POLICY_NOANONYMOUS","true");
        server = serv;
        socket = s;
        shouldClose = false;
        username = null;
        uname = null;
        haventTriedRead = 0;
        socket.setReceiveBufferSize(Server.BUFFERLEN);
        socket.setSendBufferSize(Server.BUFFERLEN);
        socket.setSoTimeout(1000000);
        instream = new BufferedInputStream(socket.getInputStream());
        outstream = socket.getOutputStream();
        outchannel = Channels.newChannel(outstream);
        bufferpos = 0;
        sasl = null;
        socket.setTcpNoDelay(true);
        taskInactivityLimit = inactivityLimit;
        if(taskInactivityLimit>0)
        	limitInactivity=true;
        lastActivity = System.currentTimeMillis();
    }
    
    /* The shouldClose method checks the last activity time (set in run method) and returns true 
    *  if past the time limit. Other methods set shouldClose after exceptions or once the task is finished. 
    */
    public boolean shouldClose() {
    	//Check if we should shut this connection down
        if(limitInactivity && System.currentTimeMillis()-lastActivity > taskInactivityLimit){
        	shouldClose=true;
        	server.getLogger().log(Level.INFO,"serverTask","connection idle more than "+taskInactivityLimit+"ms: closing connection");
        }
        
        return shouldClose;
    }
    public void close () {
        try {
            if (server.debug()) {
                System.err.println("Closing Socket " + socket + " for " + this);
            }
            if (!socket.isClosed()) {
                socket.close();
            }
        } catch (Exception e) {
            // ignore it
        }
    }
    public boolean inputAvailable() {
        boolean avail = false;
        if (bufferpos >= buffer.length) {
            shouldClose = true;
            System.err.println("inputAvailable: Buffer was full.  Closing");
            return false;
        }
        try {
            avail = instream.available() > 0;
            if (avail || haventTriedRead++ > 1000) {
                haventTriedRead = 0;
                socket.setSoTimeout(1);
                int r = instream.read();
                socket.setSoTimeout(1000000);
                if (r == -1) {
                    avail = false;
                    shouldClose = true;
                    System.err.println("inputAvailable: Connection Closed");
                } else {
                    synchronized(buffer) {
                        buffer[bufferpos++] = (byte)r;
                    }
                    //                    System.err.println("iA has read through " + new String(buffer,0,bufferpos));
                }
            }                
        } catch (SocketTimeoutException e) {
            // timeout means the socket is still open but there's no data.  That's ok.
            try {
                socket.setSoTimeout(1000000);   
            } catch (IOException e2) {
                server.getLogger().logp(Level.INFO,"serverTask","inputAvailable","socket.setSoTimeout",e2);
                avail = false;
                shouldClose = true;                
            }
            avail = false;
        } catch (IOException e) {
            e.printStackTrace();
            avail = false;
            shouldClose = true;
        }
        		
        return avail;
    }
    /** prints the response header signifying a valid request.  Only happens after
     *  the ServerTask has read enough information from the socket and done
     *  whatever else needs doing to be sure that it can satisfy the request.
     */
    public void printOK() throws IOException {
        printString("OK\n");
    }
    /** prints the response header signifying an invalid request
     */
    public void printInvalid(String reason) throws IOException {
        printString("INVALID " + reason + "\n");
    }
    /** prints the response header signifying lack of permissions */
    public void printAuthError() throws IOException {
        printString("Permission Denied\n");
    }
    /** sends the string s to the client
     */
    public void printString(String s) throws IOException {
        //         if (server.debug()) {
        //        System.err.println("SEND " + s);
        //         }
        outstream.write(s.getBytes());
        outstream.flush();
    }
    /**
     * Reads a line from the socket and returns it
     */
    public String readLine() throws IOException {
        int i = 0;
        boolean done = false;
        for (i = 0; i < bufferpos; i++) {
            if (buffer[i] == '\n') {
                String out = new String(buffer,0,i);
                synchronized(buffer) {
                    for (int j = 0; j < i; j++) {
                        buffer[j] = buffer[i+j+1];
                    }
                    bufferpos -= i+1;
                }
                //System.err.println("READ EXISTING " + out);
                return out;
            }
        }
        while (instream.available() > 0 && 
               (i = instream.read()) != -1) {
            if (i == '\n') {
                done = true;
                break;
            } else {
                if (bufferpos >= buffer.length) {
                    server.getLogger().logp(Level.WARNING,"ServerTask","readLine " + toString(),"readline would overflow.  quitting");
                    shouldClose = true;
                    return null;
                }
                buffer[bufferpos++] = (byte)i;
            }
        }
        if (i == -1) {
            shouldClose = true;
            server.getLogger().logp(Level.WARNING,"ServerTask","readLine " + toString(),"error in readline.  quitting");
            return null;
        }
        if (done) {
            String out = new String(buffer,0,bufferpos);
            bufferpos = 0;
            //System.err.println("READ " + out);
            return out;
        } else {
            //System.err.println("incomplete line in readline");
            return null;
        }
    }
    /**
     * main method for the task.  This method is asynchronous- it shouldn't block too long on the client.  It does block
     * on disk reads and such and it does block on the client while waiting, eg, for more hits to store.  It does
     * not block while reading parameters or waiting for the next command.
     *
     * The asynchronous behavior is achieved through readLine(), which returns null if it doesn't have a complete line.
     * Typically run then returns and waits until it's called again, at which point there is hopefully a complete
     * line to deal with.  aside from authenticate(), the basic procedure is
     *   - read the type of request
     *   - read the list of request parameters
     *   - call processRequest()
     *
     * Run CANNOT throw any exceptions in the current model.  If the outer block gets an exception, it swallows it
     * and returns, setting shouldClose = true to indicate that this connection to a client should be closed.
     */
    public void run() {
        try {
        	lastActivity = System.currentTimeMillis();
            if (username == null) {
                if (!authenticate()) {
                    server.getLogger().logp(Level.INFO,"serverTask","run " + toString(),"not authenticated in ");
                    printAuthError();
                    shouldClose = true;
                    return;
                }
                if (username == null) { 
                    return ;
                }
                server.getLogger().logp(Level.INFO,"ServerTask","run " + toString(), " authenticated " + username + " from " + socket.getInetAddress() + ":" + socket.getPort());
                printString("authenticated as " + username + "\n");
            }
            while (true) {
                String p = readLine();
                if (p == null) { 
                    break;
                } else {
                    if (p.equals("ENDREQUEST")) {
                        String error = request.parse(args);
                        if (error == null) {
                            processRequest();
                        } else {
                            server.getLogger().logp(Level.INFO,"ServerTask","run()" + toString(), "error parsing request: " + error);
                            printString("error parsing request: " + error + "\n");
                        }
                        args.clear();
                        if (outstream != null) { outstream.flush(); }                            
                        break;
                    } else {
                        args.add(p);
                        if (args.size() > MAXPARAMLINES) {                            
                            shouldClose = true;
                            break;
                        }
                    }
                }
            }
        } catch (Exception e) {
            server.getLogger().logp(Level.INFO,"serverTask","run " + toString(),"error " + e.toString(),e);
            args.clear();
            shouldClose = true;
            System.gc();
            System.runFinalization();
            return;
        }
    }
    /**
     * performs authentication exchange over the socket and sets the username field
     * if successful.  Returns true if authenticate should continue or is successful.
     * Returns false if authenticate has failed.
     */
    public boolean authenticate() throws IOException {
        /* The SASL client and server give you back bytes to send
           to the other side.  We achieve this by sending a length
           line first (ascii encoded integer followed by '\n')
           and then the raw bytes of the SASL exchange.  Two complexities:
           1) input.read() doesn't necessarily read the expected number
           of bytes all at once, so we have to loop around it until
           it does.
           2) I had problems with isComplete() returning true at different
           times in the client and server.  The server sends an isComplete() byte
           at the end of the loop to tell the client when it's done.
        */
        if (uname == null) {
            uname = readLine();
            if (uname == null) {
                return true;
            }
        }
        if (sasl == null) {
            sasl = Sasl.createSaslServer(Server.SaslMechanisms[0],
                                         "readdb",
                                         socket.getInetAddress().getCanonicalHostName(),
                                         saslprops,
                                         new ServerTaskCallbackHandler(server,uname));
        }
        if (sasl == null || sasl.isComplete()) {
            outstream.write("0\n".getBytes());
            outstream.write((byte)0);
            outstream.flush();
            server.getLogger().logp(Level.INFO,"ServerTask","authenticate " + toString(),"Failed Authentication for " + uname);
            return false;
        }
        while (!sasl.isComplete()) {
            try {
                String l = readLine();
                if (l == null) {
                    return true;
                }
                int length = Integer.parseInt(l);
                byte[] response = new byte[length];
                int read = 0;
                while (read < length) {
                    read += instream.read(response, read, length - read);
                }                
                byte[] challenge = sasl.evaluateResponse(response);
                if (challenge == null) {
                    challenge = new byte[0];
                }
                String s = challenge.length + "\n";
                outstream.write(s.getBytes());
                outstream.write(challenge);
                outstream.write(sasl.isComplete() ? (byte)0 : (byte)1);
                outstream.flush();
            } catch (Exception e) {
                server.getLogger().logp(Level.INFO,"serverTask","authenticate " + toString(),e.toString(),e);
                outstream.write("0\n".getBytes());
                outstream.write((byte)0);
                outstream.flush();
                break;
            }
        }
        if (sasl.isComplete() && sasl.getAuthorizationID().equals(uname)) {
            username = sasl.getAuthorizationID();
            sasl.dispose();
            return true;
        } else {
            sasl.dispose();
            server.getLogger().logp(Level.INFO,"ServerTask","authenticate " +toString(),"Failed Authentication for " + uname);
            return false;
        }

    }


    /** reads and handles a request on the Socket.
     */
    public void processRequest () {
        try {
            if (request.alignid != null) {
                Lock.readLock(request.alignid);
            }
            if (request.type.equals("ping")) {
                processPing();
            }else if (request.type.equals("exists")) {
                processExists();            
            } else if (request.type.equals("storesingle")) {
                processSingleStore();
            } else if (request.type.equals("storepaired")) {
                processPairedStore();
            } else if (request.type.equals("reindex")) {
                processReindex();
            } else if (request.type.equals("bye")) {
                shouldClose = true;
            } else if (request.type.equals("getchroms")) {
                processGetChroms();
            } else if (request.type.equals("getacl")) {
                processGetACL();            
            } else if (request.type.equals("setacl")) {
                processSetACL();
            } else if (request.type.equals("deletealign")) {
                processDeleteAlignment();
            } else if (request.type.equals("addtogroup")) {
                processAddToGroup();
            } else if (request.type.equals("shutdown")) {
                server.getLogger().logp(Level.INFO,"ServerTask","processRequest " + toString(),"Received shutdown from " + username);
                if (server.isAdmin(username)) {
                    printOK();
                    server.keepRunning(false);
                } else {
                    printAuthError();
                }
                shouldClose = true;
            } else {
                processFileRequest();
            }
        } catch (Exception e) {
            server.getLogger().logp(Level.INFO,"ServerTask","processRequest " + toString(),"Error in request " + request.toString());
            server.getLogger().logp(Level.INFO,"ServerTask","processRequest " + toString(),"Exception " + e.toString(),e);
            e.printStackTrace();
            StackTraceElement[] elts = e.getStackTrace();
            StringBuffer sb = new StringBuffer();
            for (int i = 0; i < elts.length; i++) {
                sb.append(elts[i].toString());
            }
            server.getLogger().logp(Level.INFO,"ServerTask","processRequest " + toString(),"Trace " + sb.toString());   
        } catch (AssertionError e) {
            server.getLogger().logp(Level.INFO,"ServerTask","processRequest " + toString(),"Error in request " + request.toString());
            server.getLogger().logp(Level.INFO,"ServerTask","processRequest " + toString(),"Exception " + e.toString(),e);
            e.printStackTrace();
            StackTraceElement[] elts = e.getStackTrace();
            StringBuffer sb = new StringBuffer();
            for (int i = 0; i < elts.length; i++) {
                sb.append(elts[i].toString());
            }
            server.getLogger().logp(Level.INFO,"ServerTask","processRequest " + toString(),"Trace " + sb.toString());   

        } finally {
            Lock.releaseLocks();
        }
    }
    /**
     * Handles the subset of requests that deal with a particular
     * file that we expect to exist
     */
    public void processFileRequest()  throws IOException{
        assert(request != null);
        assert(request.alignid != null);
        assert(request.chromid != null);        
        if (request.alignid == null || request.alignid.length() == 0) {
            printString("null or empty alignment " + request.alignid + "\n");
            return;
        }
        if (request.chromid == null) {
            printString("null chromosome\n");
            return;
        }
        File directory = new File(server.getAlignmentDir(request.alignid));
        if (!directory.exists()) {
            printString("No Such Alignment\n");
            return;
        } 
        AlignmentACL acl = null;        
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            // happens if the file doesn't exist or if we can't read it at the OS level
            server.getLogger().logp(Level.INFO,"ServerTask","processFileRequest "+ toString(),
                                   String.format("read error on acl for %s : %s",
                                                 request.alignid,
                                                 e.toString()));
            printInvalid(e.toString());
            return;
        }
        if (!authorizeRead(acl)) {
            server.getLogger().logp(Level.INFO,"ServerTask","processFileRequest "+toString(),
                                   String.format("%s can't read %s",
                                                 username,
                                                 request.alignid));
            printAuthError();
            return;
        }
        Header header=null;
        Hits hits=null;
        try {
            if (request.isPaired) {
                hits = server.getPairedHits(request.alignid, request.chromid, request.isLeft);
                header = server.getPairedHeader(request.alignid, request.chromid, request.isLeft);
            } else {
                hits = server.getSingleHits(request.alignid, request.chromid, request.isType2);
                header = server.getSingleHeader(request.alignid, request.chromid, request.isType2);
            }
        } catch (IOException e) {
            // We already know the alignment exists and we have permissions at this point. 
        	// This exception may happen if the file doesn't exist or if we can't read it at the OS level
        	//Let's choose to ignore the exception for now, and to have the methods below behave like there are 0 hits in this case
        	header=null;
        	hits=null;
        	
            /*server.getLogger().logp(Level.INFO,"ServerTask","processFileRequest " + toString(),
                                   String.format("read error on header or hits for %s, %d, %s : %s",
                                                 request.alignid, request.chromid, request.isLeft,
                                                 e.toString()));
            printInvalid(e.toString());
            return;
            */
        }
        if (request.type.equals("count")) {
            processCount(header,hits);
        } else if (request.type.equals("weight")) {
            processWeight(header,hits);
        } else if (request.type.equals("numpositions")) {
            processNumPositions(header,hits);
        } else if (request.type.equals("numpairpositions")) {
            processNumPairPositions(header,hits);
        } else if (request.type.equals("histogram")) {
            processHistogram(header,hits);
        } else if (request.type.equals("weighthistogram")) {
            processWeightHistogram(header,hits);
        } else if (request.type.equals("gethits")) {
            processGetHits(header,hits);
        } else if (request.type.equals("checksort")) {
            processCheckSort(header,hits);
        } else {
            printInvalid("request type");
        }
        hits = null;
        header = null;
    }

    /* returns true iff the user named in the username field is allowed to
       access this file, or if the user is in the admin group
    */
    public boolean authorizeRead(AlignmentACL acl) {
        return authorize(acl.getReadACL());
    }
    public boolean authorizeWrite(AlignmentACL acl) {
        return authorize(acl.getWriteACL());
    }
    public boolean authorizeAdmin(AlignmentACL acl) {
        return authorize(acl.getAdminACL());
    }
    private boolean authorize(Set<String> acl) {
        if (acl.contains(username)) {
            return true;
        }
        for (String g : acl) {
            if (server.groupContains(username, g)) {
                return true;
            }
        }
        if(server.isAdmin(username))
        	return true;
        return false;
    }
    public void processAddToGroup() throws IOException {
        String princ = request.map.get("princ");
        String group = request.map.get("group");
        if (!server.isAdmin(username)) {
            printAuthError();
            return;
        }
        if (princ == null || group == null) {
            printString("Must supply princ and group :" + princ + "," + group + "\n");
            return;
        }
        server.addToGroup(this,group,princ);
        printOK();
    }
    public void processGetACL() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            printString("No such alignment\n");
            return;
        }
        if (!authorizeAdmin(acl)) {
            printAuthError();
            return;
        }
        printOK();
        StringBuffer sb = new StringBuffer();
        sb.append("READ\n");
        sb.append(acl.getReadACL().size() + "\n");
        for (String s : acl.getReadACL()) {
            sb.append(s + "\n");
        }
        sb.append("WRITE\n");
        sb.append(acl.getWriteACL().size() + "\n");
        for (String s : acl.getWriteACL()) {
            sb.append(s + "\n");
        }
        sb.append("ADMIN\n");
        sb.append(acl.getAdminACL().size() + "\n");
        for (String s : acl.getAdminACL()) {
            sb.append(s + "\n");
        }
        printString(sb.toString());
    }
    public void processSetACL() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        Lock.writeLock(request.alignid);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            printString("No such alignment\n");
            return;
        }
        if (!authorizeAdmin(acl) && !server.isAdmin(username)) {
            printAuthError();
            return;
        }
        for (int i = 0; i < request.list.size(); i++) {
            String line = request.list.get(i);
            String pieces[] = line.split(" ");
            /* line format is principal [add|delete] [read|write|admin] */
            Set<String> aclset = pieces[2].equals("admin") ? acl.getAdminACL() : 
                (pieces[2].equals("write") ? acl.getWriteACL() : 
                 (pieces[2].equals("read") ? acl.getReadACL() : null));
            if (aclset == null) {
                printString("Bad ACL Type " + pieces[2] + "\n");
                continue;
            }
            if (pieces[1].equals("add")) {
                aclset.add(pieces[0]);
            } else if (pieces[1].equals("delete")) {
                aclset.remove(pieces[0]);
            } else {
                printString("Bad Operation Type " + pieces[1] + "\n");
                continue;
            }
        }
        acl.writeToFile(server.getACLFileName(request.alignid));
        server.removeACL(request.alignid);
        printString("OK\n");                        
    }
    /** Get a ping, return a pong
     */
    public void processPing() throws IOException {
        assert(request != null);
        printString("pong\n");
    }
    /** reads two lines from socket: alignment id and chromosome id.
     * returns "exists" or "unknown" to indicate whether the 
     * server knows about that pair
     */
    public void processExists() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        try {
            AlignmentACL acl = server.getACL(request.alignid);
            if (authorizeRead(acl)) {
                printString("exists\n");
            } else {
                printString("exists but no read permissions\n");
            }
        } catch (Exception e) {
            server.getLogger().logp(Level.INFO,"serverTask","processExists " +toString(),e.toString(),e);
            printString("unknown\n");
        }
    }
    /**
     * Returns the list of chromosomes for an alignment
     */
    public void processGetChroms() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            printString("No Such Alignment\n");
            return;
        }
        if (!authorizeRead(acl)) {
            printAuthError();
            return;
        }        
        Set<Integer> chroms = server.getChroms(request.alignid,
                							   request.isType2,
                                               request.isPaired,
                                               request.isLeft);
        if (chroms == null) {
            printString("No Such Alignment\n");
            return;
        }

        printOK();
        printString(chroms.size() + "\n");
        for (Integer i : chroms) {
            printString(i + "\n");
        }            
    }
    /**
     * Deletes an alignment: the header and hits files, acl file, and the directory are removed
     *
     */
    public void processDeleteAlignment() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        AlignmentACL acl = server.getACL(request.alignid);
        if (!authorizeAdmin(acl)) {
            printAuthError();
            return;
        }
        Lock.writeLock(request.alignid);
        /* step one is to de-cache all the files */
        Set<Integer> chroms = server.getChroms(request.alignid,false,request.isPaired,true);
        if (request.isPaired) {
            chroms.addAll(server.getChroms(request.alignid,false, true,false));
        }else
        	chroms.addAll(server.getChroms(request.alignid,true, false,true));
        for (int c : chroms) {
            if (request.isPaired) {
                server.removePairedHits(request.alignid, c, true);
                server.removePairedHits(request.alignid, c, false);
                server.removePairedHeader(request.alignid, c, true);
                server.removePairedHeader(request.alignid, c, false);
            } else {
                server.removeSingleHits(request.alignid, c, true);
                server.removeSingleHeader(request.alignid, c, true);
                server.removeSingleHits(request.alignid, c, false);
                server.removeSingleHeader(request.alignid, c, false);
            }
        }

        /* now do the deletes */
        boolean allDeleted = true;
        File directory = new File(server.getAlignmentDir(request.alignid));
        String prefix = directory.getCanonicalPath() + System.getProperty("file.separator");
        File[] files = directory.listFiles();
        List<String> toDelete = new ArrayList<String>();
        /* list of files to delete:
           datafiles first, then the directory itself
        */
        boolean allgone = true;
        for (int i = 0; i < files.length; i++) {
            String name = files[i].getName();
            if (name.equals("acl.txt")) { continue;}
            if (request.isPaired == null) {
                toDelete.add(prefix+name);
            } else {
                boolean pairedfile = name.indexOf(".prleft.") > 0 ||
                    name.indexOf(".prright.") > 0 ||
                    name.indexOf(".pairedleftindex") > 0 ||
                    name.indexOf(".pairedrightindex") > 0 ||
                    name.indexOf(".paircode") > 0;  
                boolean singlefile = name.indexOf("singleindex") > 0|| name.indexOf("singlet2index") > 0||
                    name.indexOf("spositions") > 0 || name.indexOf("st2positions") > 0 ||
                    name.indexOf("sweights") > 0 || name.indexOf("st2weights") > 0 ||
                    name.indexOf("slas") > 0 ||name.indexOf("st2las") > 0;
                if (request.isPaired && pairedfile) {
                    toDelete.add(prefix + name);
                } else if (!request.isPaired && singlefile) {
                    toDelete.add(prefix+name);
                } else {
                    allgone = false;
                }

            }
        }       
        if (allgone) {
            toDelete.add(server.getACLFileName(request.alignid));
            server.removeACL(request.alignid);
        }
        File f;

        for (String fname : toDelete) {
            // file system delete
            f = new File(fname);
            boolean deleted = f.delete();
            allDeleted = allDeleted && deleted;
            if (!deleted) {
                server.getLogger().logp(Level.INFO,"ServerTask","processDeleteAlignment "+ toString(),
                                       "ServerTask.processDeleteAlignment didn't delete " + fname);
            }
        }
        if(allgone){
        	// file system delete directory
            f = new File(directory.getCanonicalPath());
        	boolean deleted = f.delete();
            allDeleted = allDeleted && deleted;
            if (!deleted) {
                server.getLogger().logp(Level.INFO,"ServerTask","processDeleteAlignment "+ toString(),
                                       "ServerTask.processDeleteAlignment didn't delete directory " + directory.getCanonicalPath());
            }
        }
        
        if (allDeleted) {
            printOK();
        } else {
            printString("Partially Deleted\n");
        }
    }
    /** creates or appends to a set of hits.  
     *
     * If the chromosome file doesn't exist yet, then create a new one and dump in positions and weights.
     * If it does exist, then create a new file and merge the old file with the new
     * set of hits.
     */
    public void processSingleStore() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        assert(request.chromid != null);        
        int numHits = 0;
        try {
            numHits = Integer.parseInt(request.map.get("numhits"));
        } catch (NumberFormatException e) {
            server.getLogger().logp(Level.INFO,"ServerTask","processSingleStore "+ toString(),
                                   "Invalid numhits " + request.map.get("numhits"),e);
            printString("Invalid numhits value : " + request.map.get("numhits") + "\n");
            return;
        }
        printOK();
        if (numHits == 0) {
            printOK();
            return;
        }
        Lock.writeLock(request.alignid);

        IntBP positions = new IntBP(numHits);
        FloatBP weights = new FloatBP(numHits);
        IntBP las = new IntBP(numHits);
        ReadableByteChannel rbc = Channels.newChannel(instream);
        Bits.readBytes(positions.bb, rbc);
        Bits.readBytes(weights.bb, rbc);
        Bits.readBytes(las.bb, rbc);
        
        SingleHit[] newhits = new SingleHit[numHits];
        for (int i = 0; i < numHits; i++) {
            newhits[i] = new SingleHit(request.chromid,
                                       positions.get(i),
                                       weights.get(i),
                                       Hits.getStrandOne(las.get(i)),
                                       Hits.getLengthOne(las.get(i)));
        }
        for (int i = 1; i < newhits.length; i++) {
            if (newhits[i-1].compareTo(newhits[i]) > 0) {
                throw new RuntimeException(String.format("at %d : %d vs %d",
                                                         i, newhits[i-1].pos, newhits[i].pos));
            }
        }

        positions = null;
        weights = null;
        las = null;

        SingleHit[] hits = null;

        /* if the alignment already exists, read in the old hits */
        Set<Integer> chroms = server.getChroms(request.alignid, request.isType2, false,false);
        try {
            if (chroms != null && chroms.contains(request.chromid)) {
                try {        
                    /* sure we're allowed to write here */
                    AlignmentACL acl = server.getACL(request.alignid);
                    if (!authorizeRead(acl) || !authorizeWrite(acl)) {
                        printAuthError();
                        return;
                    }
                } catch (Exception e) {
                    server.getLogger().logp(Level.INFO,"serverTask","processSingleStore "+toString(),e.toString(),e);
                    printInvalid(e.toString());
                    return;
                }
                try {
                    server.getSingleHits(request.alignid,
                                         request.chromid, 
                                         request.isType2).appendSingleHits(newhits,
                                                                           server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                                                           request.chromid, request.isType2);
                } catch (Exception e) {
                    server.getLogger().logp(Level.INFO,"ServerTask","processSingleStore "+toString(),"error writing hits",e);
                    printInvalid(e.toString());
                    return;
                }
            } else {
                /* this is a new alignment, so set a default ACL */
                AlignmentACL acl = new AlignmentACL();
                try {
                    acl.readFromFile(server.getDefaultACLFileName());
                } catch (IOException e) {
                    // no default acl, so dont' worry.
                }
                File dir = new File(server.getAlignmentDir(request.alignid));
                if (!dir.exists() && !dir.mkdirs()) {
                    server.getLogger().logp(Level.INFO,"ServerTask","processSingleStore "+ toString(),"Can't create directories for " + request.alignid + ":" + server.getAlignmentDir(request.alignid));
                    printAuthError();
                    return;
                }
                acl.getAdminACL().add(username);
                acl.getWriteACL().add(username);
                acl.getReadACL().add(username);
                acl.writeToFile(server.getACLFileName(request.alignid));        
                server.removeACL(request.alignid); // make sure the server doesn't have this ACL cached
                SingleHits.writeSingleHits(newhits,
                                           server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                           request.chromid,
                                           request.isType2);
            }
            SingleHits singlehits = new SingleHits(server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                                   request.chromid,
                                                   request.isType2);
            Header header = new Header(singlehits.getPositionsBuffer().getib());
            header.writeIndexFile(server.getSingleHeaderFileName(request.alignid,
                                                                 request.chromid,
                                                                 request.isType2));
        } catch (IOException e) {
            server.getLogger().logp(Level.INFO,"ServerTask","processSingleStore "+ toString(),"IOException trying to save files : " + e.toString(),e);
            return;
        }
        printOK();
        server.removeSingleHits(request.alignid, request.chromid, request.isType2);
        server.removeSingleHeader(request.alignid, request.chromid, request.isType2);
    }

    public void processPairedStore() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        assert(request.chromid != null);        
        assert(request.isLeft != null);
        int numHits = 0;
        try {
            numHits = Integer.parseInt(request.map.get("numhits"));
        } catch (NumberFormatException e) {
            server.getLogger().logp(Level.INFO,"ServerTask","processPairedStore "+ toString(),"Invalid numhits " + request.map.get("numhits"),e);
            printString("Invalid numhits value : " + request.map.get("numhits") + "\n");
            return;
        }
        printOK();
        if (numHits == 0) {
            printOK();
            return;
        }

        Lock.writeLock(request.alignid);
        File f;
        f = new File(server.getAlignmentDir(request.alignid));
        if (!f.exists()) {
            if (!f.mkdirs()) {
                server.getLogger().logp(Level.INFO,"ServerTask","processPairedStore "+ toString(), "Can't create directories for " + request.alignid + ":" + server.getAlignmentDir(request.alignid));
                printAuthError();
                return;
            }
        }
        f = new File(server.getACLFileName(request.alignid));
        if (!f.exists()) {
            /* this is a new alignment, so set a default ACL */
            AlignmentACL acl = new AlignmentACL();
            try {
                acl.readFromFile(server.getDefaultACLFileName());
            } catch (IOException e) {
                // no default acl, so dont' worry.
            }
            if (!(new File(server.getAlignmentDir(request.alignid))).mkdirs()) {
            }
            acl.getAdminACL().add(username);
            acl.getWriteACL().add(username);
            acl.getReadACL().add(username);
            acl.writeToFile(server.getACLFileName(request.alignid));        
            server.removeACL(request.alignid); // make sure the server doesn't have this ACL cached
        } else {
            try {        
                /* sure we're allowed to write here */
                AlignmentACL acl = server.getACL(request.alignid);
                if (!authorizeRead(acl) || !authorizeWrite(acl)) {
                    printAuthError();
                    return;
                }
            } catch (IOException e) {
                server.getLogger().logp(Level.INFO,"ServerTask","processPairedStore "+ toString(), "IOException reading acls : " + e.toString(),e);
                printInvalid(e.toString());
                return;
            }
        }

        PairedHit[] hits = new PairedHit[numHits];
        IntBP positions = new IntBP(numHits);
        FloatBP weights = new FloatBP(numHits);
        IntBP paircodes = new IntBP(numHits);
        IntBP las = new IntBP(numHits);
        IntBP otherchrom = new IntBP(numHits);
        IntBP otherpos = new IntBP(numHits);
        ReadableByteChannel rbc = Channels.newChannel(instream);
        Bits.readBytes(positions.bb, rbc);
        Bits.readBytes(weights.bb, rbc);
        Bits.readBytes(paircodes.bb, rbc);
        Bits.readBytes(las.bb, rbc);
        Bits.readBytes(otherchrom.bb,rbc);
        Bits.readBytes(otherpos.bb,rbc);
        for (int i = 0; i < positions.limit(); i++) {
            hits[i] = new PairedHit(request.chromid,
                                    positions.get(i),
                                    Hits.getStrandOne(las.get(i)),
                                    Hits.getLengthOne(las.get(i)),
                                    otherchrom.get(i),
                                    otherpos.get(i),
                                    Hits.getStrandTwo(las.get(i)),
                                    Hits.getLengthTwo(las.get(i)),
                                    weights.get(i),
                                    paircodes.get(i));
        }
        positions = null;
        weights = null;
        paircodes=null;
        las = null;
        otherchrom = null;
        otherpos = null;

        try {
            appendPairedHits(hits,true);
            appendPairedHits(hits,false);
            hits = null;
        } catch (IOException e) {
            server.getLogger().logp(Level.INFO,"ServerTask","processPairedStore "+toString(), "IOException trying to save files : " + e.toString(),e);
            printString("Failed to write hits : " + e.toString() + "\n");
            return;
        }
        printOK();
    }

    private void appendPairedHits(PairedHit[] newhits,
                                  boolean isLeft) throws IOException {
        /* instead of the map, could have a list passed, sort, and then use subList 
           to save some memory */

        Map<Integer,List<PairedHit>> map = new HashMap<Integer,List<PairedHit>>();
        for (PairedHit h : newhits) {
            int c = isLeft ? h.leftChrom : h.rightChrom;
            if (!map.containsKey(c)) {
                map.put(c, new ArrayList<PairedHit>());
            }
            map.get(c).add(h);
        }
        newhits = null;
        Comparator<PairedHit> comp = isLeft ? new PairedHitLeftComparator() : new PairedHitRightComparator();
        for (int chromid : map.keySet()) {
            List<PairedHit> nhlist = map.get(chromid);
            Collections.sort(nhlist, comp);
            try {
                PairedHits oldhits = server.getPairedHits(request.alignid,
                                                          chromid,
                                                          isLeft);
                oldhits.appendPairedHits(nhlist,server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"), chromid, isLeft);
            } catch (FileNotFoundException e) {
                PairedHits.writePairedHits(nhlist, server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"), chromid, isLeft);
            } 
            PairedHits pairedhits = new PairedHits(server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                                   chromid, 
                                                   isLeft);
            Header header = new Header(pairedhits.getPositionsBuffer().getib());
            header.writeIndexFile(server.getPairedHeaderFileName(request.alignid,
                                                                 chromid,
                                                                 isLeft));
            server.removePairedHits(request.alignid, chromid, isLeft);
            server.removePairedHeader(request.alignid, chromid, isLeft);
        }
    }

    public void processReindex() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        assert(request.chromid != null);        
        Lock.writeLock(request.alignid);
        if (request.isPaired) {
            PairedHits hits = server.getPairedHits(request.alignid, request.chromid, true);
            Header header = new Header(hits.getPositionsBuffer().getib());
            header.writeIndexFile(server.getPairedHeaderFileName(request.alignid,
                                                                 request.chromid,
                                                                 true));            
            server.removePairedHeader(request.alignid, request.chromid,true);

            hits = server.getPairedHits(request.alignid, request.chromid, false);
            header = new Header(hits.getPositionsBuffer().getib());
            header.writeIndexFile(server.getPairedHeaderFileName(request.alignid,
                                                                 request.chromid,
                                                                 false));            
            server.removePairedHeader(request.alignid, request.chromid,false);

        } else {
            SingleHits hits = server.getSingleHits(request.alignid, request.chromid, request.isType2);
            Header header = new Header(hits.getPositionsBuffer().getib());
            header.writeIndexFile(server.getSingleHeaderFileName(request.alignid,
                                                                 request.chromid,
                                                                 request.isType2));
            server.removeSingleHeader(request.alignid, request.chromid, request.isType2);       
        }
        printOK();
    }

    public void processCount(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
    		printOK();
    		printString("0\n");
    	}else{
	        printOK();
	        if (request.start == null && request.end == null && request.minWeight == null && request.isPlusStrand == null) {
	            printString(Integer.toString(header.getNumHits()) + "\n");
	            return;
	        }
	        if (request.start == null) {
	            request.start = 0;
	        }
	        if (request.end == null) {
	            request.end = Integer.MAX_VALUE;
	        }
	        int first = header.getFirstIndex(request.start == null ? 0 : request.start);
	        int last = header.getLastIndex(request.end == null ? Integer.MAX_VALUE : request.end);
	        printString(Integer.toString(hits.getCountBetween(first,last,request.start,request.end,request.minWeight, request.isPlusStrand)) + "\n");
    	}
    }
    public void processWeight(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
    		printOK();
    		printString("0.0\n");
    	}else{
    		printOK();
	        if (request.start == null) {
	            request.start = 0;
	        }
	        if (request.end == null) {
	            request.end = Integer.MAX_VALUE;
	        }
	        int first = header.getFirstIndex(request.start);
	        int last = header.getLastIndex(request.end);
	        printString(Double.toString(hits.getWeightBetween(first,last,request.start,request.end,request.minWeight, request.isPlusStrand)) + "\n");
    	}
    }
    public void processNumPositions(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
    		printOK();
    		printString("0\n");
    	}else{
	        printOK();
	        if (request.start == null) {
	            request.start = 0;
	        }
	        if (request.end == null) {
	            request.end = Integer.MAX_VALUE;
	        }
	        int first = header.getFirstIndex(request.start == null ? 0 : request.start);
	        int last = header.getLastIndex(request.end == null ? Integer.MAX_VALUE : request.end);
	        printString(Integer.toString(hits.getNumPositionsBetween(first,last,request.start,request.end,request.minWeight, request.isPlusStrand)) + "\n");
    	}
    }
    public void processNumPairPositions(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
    		printOK();
    		printString("0\n");
    	}else{
	        printOK();
	        if (request.start == null) {
	            request.start = 0;
	        }
	        if (request.end == null) {
	            request.end = Integer.MAX_VALUE;
	        }
	        int first = header.getFirstIndex(request.start == null ? 0 : request.start);
	        int last = header.getLastIndex(request.end == null ? Integer.MAX_VALUE : request.end);
	        if (request.isPaired)
	        	printString(Integer.toString(((PairedHits)hits).getNumPairedPositionsBetween(first,last,request.start,request.end,request.minWeight)) + "\n");
	        else
	        	printString("0\n");
    	}
    }
    public void processGetHits(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
    		printOK();
    		printString("0\n");
    	}else{
	    	int count;
	        if (!(hits instanceof PairedHits)) {
	            if (request.map.containsKey("wantotherchroms") ||
	                request.map.containsKey("wantotherpositions") ||
	                request.map.containsKey("wantpaircodes")) {
	                printString("invalid columns requested for single-ended data\n");
	                return;
	            }
	        }
	        if (request.start == null) {
	            request.start = 0;
	        }
	        if (request.end == null) {
	            request.end = Integer.MAX_VALUE;
	        }
	        int first = header.getFirstIndex(request.start);
	        int last = header.getLastIndex(request.end);
	        if (request.start == 0 && request.end == Integer.MAX_VALUE && request.minWeight == null && request.isPlusStrand == null) {
	            count = header.getNumHits();
	        } else {
	            count = hits.getCountBetween(first,last,request.start,request.end,request.minWeight, request.isPlusStrand);
	        }
	        printOK();
	        printString(Integer.toString(count) + "\n");
	        if (request.map.containsKey("wantpositions")) {
	            IntBP p = hits.getHitsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
	            Bits.sendBytes(p.bb, outchannel);
	        }
	        if (request.map.containsKey("wantweights")) {
	            FloatBP p = hits.getWeightsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
	            Bits.sendBytes(p.bb, outchannel);
	        }
	        if (request.map.containsKey("wantpaircodes")) {
	        	IntBP p = ((PairedHits)hits).getPairCodesBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
	            Bits.sendBytes(p.bb, outchannel);
	        }
	        if (request.map.containsKey("wantlengthsandstrands")) {
	            IntBP p = hits.getLASBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
	            Bits.sendBytes(p.bb, outchannel);
	        }
	        if (request.map.containsKey("wantotherchroms")) {
	            IntBP p = ((PairedHits)hits).getOtherChromsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
	            Bits.sendBytes(p.bb, outchannel);
	        }
	        if (request.map.containsKey("wantotherpositions")) {
	            IntBP p = ((PairedHits)hits).getOtherPositionsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
	            Bits.sendBytes(p.bb, outchannel);
	        }
    	}
    }
    public void processHistogram(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
    		printOK();
    		printString("0\n");
    	}else{
    		int binsize = 10;
	    	if (request.start == null) {
	            IntBP ib = hits.getPositionsBuffer();
	            request.start = ib.get(0);
	        }
	        if (request.end == null) {
	            IntBP ib = hits.getPositionsBuffer();
	            request.end = ib.get(ib.limit()-1);
	        }
	        try {
	            binsize = Integer.parseInt(request.map.get("binsize"));
	        } catch (Exception e) {
	            server.getLogger().logp(Level.INFO,"ServerTask","processHistogram "+toString(), "Exception parsing binsize : " + request.map.get("binsize"),e);
	            printString("missing or invalid bin size : " + request.map.get("binsize") + "\n");
	            return;
	        }
	        int dedup = 0, extension=0;
	        if (request.map.containsKey("dedup")) {
	            dedup = Integer.parseInt(request.map.get("dedup"));
	        }
	        if(request.map.containsKey("extension")) {
	        	extension = Integer.parseInt(request.map.get("extension"));
	        }
	        int first = header.getFirstIndex(request.start);
	        int last = header.getLastIndex(request.end);
	        int[] raw = hits.histogram(first,
	                                   last,
	                                   request.start,
	                                   request.end,
	                                   binsize,
	                                   dedup,
	                                   request.minWeight,
	                                   request.isPlusStrand,
	                                   extension);
	        int n = 0;
	        for (int i = 0; i< raw.length; i++) {
	            if (raw[i] > 0) {
	                n++;
	            }
	        }
	        int[] hist = new int[n*2];
	        int pos = 0;
	        for (int i = 0; i< raw.length; i++) {
	            if (raw[i] > 0) {
	                hist[pos*2] = request.start + binsize * i + binsize / 2;
	                hist[pos*2+1] = raw[i];
	                pos++;
	            }
	        }
	        printOK();
	        printString(Integer.toString(hist.length) + "\n");
	        Bits.sendInts(hist, outstream, buffer);
    	}
    }

    /* returns a histogram of hit weights in a region.  Inputs
     * startposition, stopposition, binsize.  Each bin's
     * value is the total amount of weight in the bin.
     *
     * bins with zero count are not included.
     */
    public void processWeightHistogram(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
    		printOK();
    		printString("0\n");
    	}else{
	    	int binsize = 10;
	        if (request.start == null) {
	            IntBP ib = hits.getPositionsBuffer();
	            request.start = ib.get(0);
	        }
	        if (request.end == null) {
	            IntBP ib = hits.getPositionsBuffer();
	            request.end = ib.get(ib.limit()-1);
	        }
	        try {
	            binsize = Integer.parseInt(request.map.get("binsize"));
	        } catch (Exception e) {
	            server.getLogger().logp(Level.INFO,"ServerTask","processWeightHistogram "+toString(), "Exception parsing binsize : " + request.map.get("binsize"),e);
	            printString("missing or invalid bin size : " + request.map.get("binsize") + "\n");
	            return;
	        }
	        int dedup = 0, extension=0;
	        if (request.map.containsKey("dedup")) {
	            dedup = Integer.parseInt(request.map.get("dedup"));
	        }
	        if(request.map.containsKey("extension")) {
	        	extension = Integer.parseInt(request.map.get("extension"));
	        }
	        int first = header.getFirstIndex(request.start);
	        int last = header.getLastIndex(request.end);
	        float[] raw = hits.weightHistogram(first,
	                                           last,
	                                           request.start,
	                                           request.end,
	                                           binsize,
	                                           dedup,
	                                           request.minWeight,
	                                           request.isPlusStrand,
	                                           extension);
	        int n = 0;
	        for (int i = 0; i< raw.length; i++) {
	            if (raw[i] > 0) {
	                n++;
	            }
	        }
	        int[] parray = new int[n];
	        float[] farray = new float[n];
	        int pos = 0;
	        for (int i = 0; i< raw.length; i++) {
	            if (raw[i] > 0) {
	                parray[pos] = request.start + binsize * i + binsize / 2;
	                farray[pos] = raw[i];
	                pos++;
	            }
	        }
	        printOK();
	        printString(Integer.toString(parray.length) + "\n");
	        Bits.sendInts(parray, outstream, buffer);        
	        Bits.sendFloats(farray, outstream, buffer);
    	}
    }
    public void processCheckSort(Header header, Hits hits) throws IOException {
    	if(header==null || hits==null){
            printString("File does not exist for this chromosome");
    	}else{
	    	IntBP ints = hits.getPositionsBuffer();
	        boolean needsort = false;
	        for (int i = 1; i < ints.limit(); i++) {
	            if (ints.get(i-1) > ints.get(i)) {
	                //                printString(String.format("Bad sort at %d : %d > %d.\n",
	                //                                          i,ints.get(i-1),ints.get(i)));
	                needsort = true;
	            }
	        }
	        if (needsort) {
	            if (hits instanceof SingleHits) {
	                server.getLogger().logp(Level.INFO,"ServerTask","processCheckSort",String.format("Resorting %s %d",request.alignid, request.chromid));
	                ((SingleHits)hits).resort(server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
	                                          request.chromid, 
	                                          request.isType2);
	                
	                server.removeSingleHits(request.alignid, request.chromid, request.isType2);
	                server.removeSingleHeader(request.alignid, request.chromid, request.isType2);       
	                hits = server.getSingleHits(request.alignid, request.chromid, request.isType2);
	                
	                header = new Header(hits.getPositionsBuffer().getib());
	                header.writeIndexFile(server.getSingleHeaderFileName(request.alignid,
	                                                                     request.chromid,
	                                                                     request.isType2));
	
	            } else {
	                printString("Can't resort paired hits");
	                return;
	            }
	        }
	        printOK();
    	}
    }
    public String toString() {
        return String.format("thread %s, user %s, remote %s:%d",
                             Thread.currentThread().toString(),
                             username,
                             socket.getInetAddress(), socket.getPort());
    }
}

/** 
 * SASL callback handler for the authenticate() method.  Provides
 * a username and server name to the CallBack
 */
class ServerTaskCallbackHandler implements CallbackHandler {
    private String uname;
    private Server server;
    public ServerTaskCallbackHandler(Server s, String u) {
        uname = u;
        server = s;
    }

    public void handle(Callback[] callbacks) {
        for (int i = 0; i < callbacks.length; i++) {
            if (callbacks[i] instanceof PasswordCallback) {
                PasswordCallback pc = (PasswordCallback)callbacks[i];
                try {
                    pc.setPassword(server.getPassword(uname).toCharArray());
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (NullPointerException e) {
                    System.err.println("No password for " + uname);
                }

            }
            if (callbacks[i] instanceof NameCallback) {
                NameCallback nc = (NameCallback)callbacks[i];
                nc.setName(uname);
            }            
            if (callbacks[i] instanceof AuthorizeCallback) {
                AuthorizeCallback ac = (AuthorizeCallback)callbacks[i];
                ac.setAuthorized(ac.getAuthenticationID().equals(ac.getAuthorizationID()));
            }

        }
    }
}
