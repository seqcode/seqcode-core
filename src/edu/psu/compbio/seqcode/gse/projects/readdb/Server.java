package edu.psu.compbio.seqcode.gse.projects.readdb;

import java.net.*;
import java.util.*;
import java.util.logging.*;
import java.io.*;
import org.apache.commons.cli.*;

/**
 * <p>ReadDB server class.  Provides configuration-type information, but
 * pushes all the request handling off onto ServerTask. Usage:
 * <pre>java edu.psu.compbio.seqcode.gse.projects.readdb.Server --datadir /path/to/data</pre>
 *
 * <p>Optional parameters
 * <ul>
 * <li>--port 52000     port to listen on
 * <li>--threads 5      number of threads to start to handle client requests
 * <li>--cachesize 100  number of chromosomes to keep files open for
 * <li>--maxconn 1000    maximum number of client connections
 * <li>--sleepiness 2   how sleepy the server should be waiting for input.  Lower values use more CPU but improve responsiveness
 * <li>--idlelimit 24  number of hours after which idle task connections are closed 
 * <li>--help           print the usage message and exit
 *
 */
public class Server {


    public static final String SaslMechanisms[] = {"CRAM-MD5","DIGEST-MD5"};

	private Logger logger;
    private int port;
    private int numThreads, cacheSize, maxConnections, sleepiness, taskIdleLimit;
    private boolean debug;
    /* topdir is the top-level directory for our data files.
      pwfile is "${topdir}/users.txt" and groupfile is 
      "${topdir}/groups.txt"
    */
    private String topdir, pwfile, groupfile;
    private boolean keepRunning;
    private Dispatch dispatch;
    private Map<String,Set<String>> groups;
    // BUFFERLEN should be a multiple of 8 to avoid problems with partial ints, floats, or doubles
    // in buffers when the buffer is allocated in bytes.
    public static final int BUFFERLEN = 8192 * 16;

    private LRUCache<Header> singleHeaders;
    private LRUCache<Header> pairedHeaders;
    private LRUCache<SingleHits> singleHits;
    private LRUCache<PairedHits> pairedHits;
    private LRUCache<AlignmentACL> acls;    

    private ServerSocket socket;

    public Server () {
        port = 52000;
        sleepiness = 4;
        numThreads = 5;
        cacheSize = numThreads * 20;
        maxConnections = 1000;
        taskIdleLimit = 24;
        topdir = "/tmp";
        keepRunning = true;
        logger = Logger.getLogger("edu.psu.compbio.seqcode.gse.projects.readdb.Server");
        logger.log(Level.INFO,"created Server");        

    }
    public void parseArgs(String[] args) throws ParseException {
        Options options = new Options();
        options.addOption("p","port",true,"port to listen on");
        options.addOption("t","threads",true,"number of threads to spawn");
        options.addOption("d","datadir",true,"directory to use for data");
        options.addOption("D","debug",false,"provide debugging output");
        options.addOption("C","cachesize",true,"how many files to keep open (this value times up to 18)");
        options.addOption("M","maxconn",true,"how many connections are allowed");
        options.addOption("S","sleepiness",true,"how sleepy the server should be while waiting for input.  1-100");
        options.addOption("L","idlelimit",true,"number of hours after which to close idle connections. Negative sets no limit.");
        options.addOption("h","help",false,"print help message");
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse( options, args, false );            
        if (line.hasOption("help")) {
            printHelp();
            System.exit(1);
        }

        if (line.hasOption("port")) {
            port = Integer.parseInt(line.getOptionValue("port"));
        }
        if (line.hasOption("threads")) {
            numThreads = Integer.parseInt(line.getOptionValue("threads"));
            cacheSize = 10 * numThreads;
        }
        if (line.hasOption("datadir")) {
            topdir = line.getOptionValue("datadir");
        }
        if (line.hasOption("cachesize")) {
            cacheSize = Integer.parseInt(line.getOptionValue("cachesize"));
        }
        if (line.hasOption("maxconn")) {
            maxConnections = Integer.parseInt(line.getOptionValue("maxconn"));
        }
        if (line.hasOption("sleepiness")) {
            sleepiness = Integer.parseInt(line.getOptionValue("sleepiness"));
            if (sleepiness < 1) {
                sleepiness = 1;
            }
            if (sleepiness > 100) {
                sleepiness = 100;
            }
        }
        if (line.hasOption("idlelimit")) {
            taskIdleLimit = Integer.parseInt(line.getOptionValue("idlelimit"));
        }

        singleHits = new LRUCache<SingleHits>(cacheSize);
        pairedHits = new LRUCache<PairedHits>(cacheSize);
        singleHeaders = new LRUCache<Header>(cacheSize);
        pairedHeaders = new LRUCache<Header>(cacheSize);
        acls = new LRUCache<AlignmentACL>(cacheSize);
        debug = line.hasOption("debug");
        logger.log(Level.INFO,String.format("Server parsed args: port %d, threads %d, directory %s",port,numThreads,topdir));
        pwfile = topdir + System.getProperty("file.separator") + "users.txt";
        groupfile = topdir + System.getProperty("file.separator") + "groups.txt";
    }
    public void printHelp() {
        System.out.println("ReadDB server process");
        System.out.println("usage: java edu.psu.compbio.seqcode.gse.projects.readdb.Server --datadir /path/to/datadir --port 52000");
        System.out.println(" [--threads 5]   use this number of worker threads to process requests.");
        System.out.println(" [--cachesize 100]  number of datasets to keep open.  Actual number of open files will be");
        System.out.println("                  up to 18 times this value");
        System.out.println(" [--maxconn 250]   maximum number of open connections");
        System.out.println(" [--debug]  print debugging output");
        System.out.println(" [--sleepiness 4]  (1-100) higher values use less CPU when idle but may incur more delay in processing requests");
    }
    public static void main(String args[]) throws Exception {
        Server server = new Server();
        server.parseArgs(args);
        server.readAndProcessGroupsFile();
        server.listen();
        System.exit(0);
    }
    public boolean keepRunning() {
        return keepRunning;
    }
    public void keepRunning(boolean k) {
        keepRunning = k;
        if (keepRunning == false && socket != null) {
            try {
                socket.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }
    public boolean debug() {return debug;}
    public int getSleepiness() {return sleepiness;}
    public void listen() throws IOException {
        Thread t = new Thread(new CacheGCHook(logger));
        t.start();
        dispatch = new Dispatch(this,numThreads, maxConnections);
        t = new Thread(dispatch);
        t.start();
        socket = new ServerSocket(port);
        socket.setReuseAddress(true);
        socket.setReceiveBufferSize(BUFFERLEN);
        socket.setSoTimeout(1000*3600*24);
        while (keepRunning) {
            try {
                Socket s = socket.accept();
                logger.log(Level.INFO,"accepted from " + s.getInetAddress());
                s.setSoLinger(false,0);
                ServerTask st = new ServerTask(this,s, taskIdleLimit*1000*3600);
                if (debug) {
                    System.err.println("New Task is " + st);
                }
                dispatch.addWork(st);
            }catch (SocketTimeoutException ste){
            	System.err.println("Connection timeout: "+ste.getMessage());
            }catch (SocketException se){
            	System.err.println("Connection error: "+se.getMessage());
        	}catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    public Logger getLogger() {return logger;}
    public String getTopDir() {
        return topdir; 
    }
    public String cleanStringForFilename(String i) {
        return i.replaceAll("[^A-Za-z0-9_\\-\\+]","_");
    }
    /* where to find the header file for this alignment
     */
    public String getAlignmentDir(String alignID) {
        alignID = cleanStringForFilename(alignID);
        return getTopDir() + System.getProperty("file.separator") + 
            alignID;
    }
    public String getACLFileName(String alignID) {
        return getAlignmentDir(alignID) + System.getProperty("file.separator") + "acl.txt";
    }
    public String getDefaultACLFileName() {
        return getTopDir() + "defaultACL.txt";
    }    
    public String getSingleHeaderFileName(String alignID,
                                          int chromID,
                                          boolean isType2) {
        return getAlignmentDir(alignID) + System.getProperty("file.separator") + chromID + (isType2 ? ".singlet2index" : ".singleindex");
    }    
    public String getPairedHeaderFileName(String alignID,
                                          int chromID,
                                          boolean isLeft) {
        return getAlignmentDir(alignID) + System.getProperty("file.separator") + chromID + ".paired" + (isLeft ? "left" : "right") + "index";
    }
    public Set<Integer> getChroms(String alignID,
                                 boolean isType2,
    							 boolean isPaired,
                                 boolean isLeft) {
        File directory = new File(getAlignmentDir(alignID));
        File[] files = directory.listFiles();
        if (files == null) {
            return null;
        }
        Set<Integer> output = new HashSet<Integer>();
        String suffix = isPaired ? (isLeft ? ".pairedleftindex" : ".pairedrightindex") : (isType2 ? ".singlet2index" : ".singleindex");
        for (int i = 0; i < files.length; i++) {
            int index = files[i].getName().indexOf(suffix);
            if (index > 0) {
                output.add(Integer.parseInt(files[i].getName().substring(0,index)));
            }
        }
        return output;
    }

    /**
     * Returns the requested Hits object.  Creates it or retrieves from cache.
     * Client code is responsible for locking the file as necessary.
     */
    public SingleHits getSingleHits(String alignID,
                                    int chrom,
                                    boolean isType2) throws IOException, SecurityException, FileNotFoundException {
        String key = alignID + chrom + isType2;
        SingleHits output = singleHits.get(key);
        if (output == null) {
            String prefix = getAlignmentDir(alignID) + System.getProperty("file.separator");
            output = new SingleHits(prefix,chrom,isType2);
            singleHits.add(key, output);
        }
        return output;
    }
    public PairedHits getPairedHits(String alignID,
                                    int chrom,
                                    boolean isLeft) throws IOException, SecurityException, FileNotFoundException {
        String key = alignID + chrom + isLeft;
        PairedHits output = pairedHits.get(key);
        if (output == null) {
            String prefix = getAlignmentDir(alignID) + System.getProperty("file.separator");
            output = new PairedHits(prefix, chrom, isLeft);
            pairedHits.add(key, output);
        }
        return output;
    }
    /**
     * Returns the requested Header object.  Creates it or retrieves from cache.
     * Client code is responsible for locking the file as necessary.
     */
    public Header getSingleHeader(String alignID, int chromID, boolean isType2) throws IOException {
        String key = alignID + chromID +isType2;
        Header output = singleHeaders.get(key);
        if (output == null) {
            output = Header.readIndexFile(getSingleHeaderFileName(alignID,chromID, isType2));
            singleHeaders.add(key, output);
        }
        return output;
    }
    public Header getPairedHeader(String alignID, int chromID, boolean isLeft) throws IOException {
        String key = alignID + chromID + isLeft;
        Header output = pairedHeaders.get(key);
        if (output == null) {
            output = Header.readIndexFile(getPairedHeaderFileName(alignID,chromID,isLeft));
            pairedHeaders.add(key, output);
        }
        return output;
    }
    /**
     * Returns the requested ACL object.  Creates it or retrieves from cache.
     * Client code is responsible for locking the file as necessary.
     */
    public AlignmentACL getACL(String alignID) throws IOException {
        AlignmentACL output = acls.get(alignID);
        if (output == null) {
            output = new AlignmentACL(getACLFileName(alignID));
            acls.add(alignID,output);
        }
        return output;
    }
    public void removeSingleHits(String alignID, int chromID, boolean isType2) {
        singleHits.remove(alignID + chromID+ isType2);
    }
    public void removePairedHits(String alignID, int chromID, boolean isLeft) {
        pairedHits.remove(alignID + chromID + isLeft);
    }
    public void removeSingleHeader(String alignID, int chromID, boolean isType2) {
        singleHeaders.remove(alignID + chromID + isType2);
    }
    public void removePairedHeader(String alignID, int chromID, boolean isLeft) {
        pairedHeaders.remove(alignID + chromID + isLeft);
    }
    public void removeACL(String alignID) {acls.remove(alignID);}
    protected void printCacheContents() {
        singleHeaders.printKeys();
        pairedHeaders.printKeys();
        singleHits.printKeys();
        pairedHits.printKeys();
        acls.printKeys();
    }

    /**
     * Returns true iff this princ is allowed
     * to create alignments.
     *
     * Currently implemented as members of the 
     * cancreate group
     */
    public boolean canCreate(String username) {
        return groupContains(username, "cancreate");
    }
    /**
     * Returns true iff this princ is a server admin.
     * Server admins can shut the server down and have read/write/admin access to any data.
     *
     * Currently implemented as members of the admin
     * group.
     */
    public boolean isAdmin(String username) {
        return groupContains(username, "admin");
    }

    /* Authenticate stuff for reading group and password files */

    public boolean groupContains(String username, String group) {
        if (!groups.containsKey(group)) {
            return false;
        }
        return groups.get(group).contains(username);
    }
    /**
     * The groups file has lines of the form
     *  groupname:user1 user2 user3 @othergroup user4
     *
     * Group names that are also usernames are removed 
     * since ACLs can contain either users or groups and we don't
     * want confusion.
     *
     */
    public Map<String,Set<String>> readGroupsFile() throws IOException {
        Map<String,Set<String>> output = new HashMap<String,Set<String>>();
        File f = new File(groupfile);
        if (!f.exists()) {
            logger.log(Level.WARNING,"No Groups file found: " + groupfile);
            throw new IOException("No Groups file found");
        }
        BufferedReader reader = new BufferedReader(new FileReader(groupfile));
        String line = null;
        while ((line = reader.readLine()) != null ) {
            String pieces[] = line.split("\\s*\\:\\s*");
            String groupname = pieces[0];
            if (pieces.length == 1) {  continue;}  // empty group
            pieces = pieces[1].split("\\s+");
            Set<String> members = output.containsKey(groupname) ? output.get(groupname) : new HashSet<String>();
            for (int i = 0; i < pieces.length; i++) {
                if (!members.contains(pieces[i])) {
                    members.add(pieces[i]);
                }
            }
            output.put(groupname,members);
        }
        reader.close();
        return output;
    }
    public void readAndProcessGroupsFile() throws IOException {
        groups = readGroupsFile();
        expand(groups, new HashMap<String,Set<String>>());
        for (String u : getUserNames()) {
            groups.remove(u);
        }
    }
    public void addToGroup(ServerTask t, String group, String princ) throws IOException {
        System.err.println("Adding " + princ + " to " + group);
        Lock.readLock("group.txt");
        Map<String,Set<String>> rawgroups = readGroupsFile();
        if (!rawgroups.containsKey(group)) {
            Set<String> s = new HashSet<String>();
            s.add(princ);
            rawgroups.put(group, s);
        }        
        if (!rawgroups.get(group).contains(princ)) {
            rawgroups.get(group).add(princ);
        }        
        System.err.println("got write lock");
        File gfile = File.createTempFile("tmp",".groups");
        PrintWriter pw = new PrintWriter(gfile);
        for (String g : rawgroups.keySet()) {
            pw.print(g + ":");
            for (String p : rawgroups.get(g)) {
                pw.print(" " + p);
            }
            pw.println();
        }
        System.err.println("Done.  rereading");
        pw.close();
        Lock.writeLock("group.txt");
        gfile.renameTo(new File(groupfile));
        readAndProcessGroupsFile();
    }
    /**
     * processes recursive group memberships
     */
    private void expand(Map<String, Set<String>> toexpand,
                        Map<String, Set<String>> expanded) {
        boolean keepgoing = false;
        for (String g : toexpand.keySet()) {
            if (!expanded.containsKey(g)) {
                expanded.put(g,new HashSet<String>());
            }
            for (String m : toexpand.get(g)) {
                if (m.matches("^@")) {
                    toexpand.get(g).remove(m);
                    String othergroup = m.replaceAll("^@","");
                    if (!expanded.get(g).contains(m)) {
                        expanded.get(g).add(m);
                        if (!toexpand.get(g).containsAll(toexpand.get(othergroup))) {
                            toexpand.get(g).addAll(toexpand.get(othergroup));                           
                        }
                    }
                }
            }
            for (String m : toexpand.get(g)) {
                if (m.matches("^@") && !expanded.get(g).contains(m)) {
                    keepgoing = true;
                }
            }
        }
        if (keepgoing) {
            logger.log(Level.INFO,"Recursing in expand");
            expand(toexpand,expanded);
        }
    }

    /** returns the password for the specified user, or 
     * null if the user is unknown
     *
     * The users file has lines of the form username:password
     */
    protected String getPassword(String username) throws IOException {
        File f = new File(pwfile);
        if (!f.exists()) {
            logger.log(Level.SEVERE,"Don't see password file " + pwfile);
        }

        BufferedReader reader = new BufferedReader(new FileReader(pwfile));
        String line = null;
        while ((line = reader.readLine()) != null ) {
            String pieces[] = line.split("\\s*\\:\\s*");
            if (pieces[0].equals(username)) {
                reader.close();
                return pieces[1];
            }
        }
        reader.close();
        return null;
    }
    /** returns the set of all known usernames
     */
    protected Set<String> getUserNames() throws IOException {
        HashSet<String> output = new HashSet<String>();
        BufferedReader reader = new BufferedReader(new FileReader(pwfile));
        String line = null;
        while ((line = reader.readLine()) != null ) {
            String pieces[] = line.split("\\s*\\:\\s*");
            output.add(pieces[0]);
        }
        return output;
    }

        

}

