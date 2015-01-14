package edu.psu.compbio.seqcode.gse.utils.database;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.tomcat.jdbc.pool.DataSource;
import org.apache.tomcat.jdbc.pool.PoolProperties;

/**
 * <code>DatabaseConnectionManager</code> is the class that all code should use to obtain database connections. 
 * Provides database connections for *roles* and uses Tomcat JDBC to manage a connection pool.
 * A role is a database resource that you want to access, eg, <b>ucsc_mm9</b> for 
 * annotations for mm9 or <b>core</b> for the core db.
 *
 * @author Shaun Mahony
 * @version 1.0
 */
public abstract class DatabaseConnectionManager {

	static Pattern oraclePattern = Pattern.compile(".*oracle.*",Pattern.CASE_INSENSITIVE);
    static Pattern mysqlPattern = Pattern.compile(".*mysql.*",Pattern.CASE_INSENSITIVE);
    public static int ORACLE = 1, MYSQL = 2, POSTGRES = 3, SQLLITE = 4, UNKNOWN;

    /*
     * Static variables
     */
    private static List<String> usedRoles = new ArrayList<String>();
    private static Map<String,String>defaultUsers = new HashMap<String,String>();
    private static Map<String,String>defaultSchemas = new HashMap<String,String>();
    private static Map<Connection, DataSource> cxnSource = new HashMap<Connection,DataSource>();
    private static Map<String,DataSource> sources = new HashMap<String, DataSource>();
    private static Map<DataSource, Integer> sourceType = new HashMap<DataSource, Integer>();
    
    /**
     * returns the default read-only connection for this role 
     */
    public static Connection getConnection(String role) throws SQLException, UnknownRoleException {
        Properties props = null;
        try {
        	role = getRealRole(role);
        	if(!usedRoles.contains(role))
        		usedRoles.add(role);
            String user, schema;
            if (defaultUsers.containsKey(role) && 
                defaultSchemas.containsKey(role)) {
                user = defaultUsers.get(role);
                schema = defaultSchemas.get(role);
            } else {
                props = getPropertiesForRole(role);
                user = props.getProperty("user");
                schema = props.getProperty("schema");
                defaultUsers.put(role,user);
                defaultSchemas.put(role,schema);
            }
            String key = role + user + schema;
            if (!sources.containsKey(key)) {
                if (props == null)
                    props = getPropertiesForRole(role);
                addDataSource(key,props);
            }
            DataSource src = sources.get(key);            
            Connection cxn = src.getConnection();
            cxnSource.put(cxn,src);
            return cxn;
        } catch (IOException ex) {
            throw new RuntimeException("Couldn't read properties for " + role,ex);
        }
    }

    /**
     * Return the username associated with the role
     * @param role
     * @return
     */
    public static String getUsername(String role){
    	Properties props = null;
        try {
            role = getRealRole(role);
            String user;
            if (defaultUsers.containsKey(role) ) {
                user = defaultUsers.get(role);
            } else {
                props = getPropertiesForRole(role);
                user = props.getProperty("user");
                defaultUsers.put(role,user);
            }
            return user;
        } catch (IOException ex) {
            throw new RuntimeException("Couldn't read properties for " + role,ex);
        }
    }
    
    /** for a given role name, eg chipchip, look up
     * CHIPCHIPROLE in the environment.  If it's set,
     * use its value as the return value which is the
     * actual role to use.
     *
     * For example, you might set CHIPCHIPROLE=finkchipchip
     * to connect to a chipchip schema that isn't the one
     * in your .chipchip_passwd
     */
    public static String getRealRole(String role) {
        String envkey = role.toUpperCase() + "ROLE";
        String val = System.getenv(envkey);
        if (val == null) {
            return role;
        } else {
            return val;
        }
    }

    public static void addDataSource(String key, Properties props) throws SQLException{
        String cs = props.getProperty("jdbcconnectstring");
        Matcher o = oraclePattern.matcher(cs);
        Matcher m = mysqlPattern.matcher(cs);
        int currType;
        
        PoolProperties p = new PoolProperties();
        p.setUrl(props.getProperty("jdbcconnectstring"));
        if (o.matches()) {
        	p.setDriverClassName("oracle.jdbc.OracleDriver");
        	currType=ORACLE;
        } else if (m.matches()) {
        	p.setDriverClassName("com.mysql.jdbc.Driver");
        	currType=MYSQL;
        } else {
            throw new SQLException("Unknown database type in " + cs);
        }
        p.setUsername(props.getProperty("user"));
        p.setPassword(props.getProperty("passwd"));
        p.setJmxEnabled(true);
        p.setTestWhileIdle(false);
        p.setTestOnBorrow(true);
        p.setValidationQuery("SELECT 1");
        p.setTestOnReturn(false);
        p.setValidationInterval(30000);
        p.setTimeBetweenEvictionRunsMillis(30000);
        p.setMaxActive(100);
        p.setInitialSize(10);
        p.setMaxWait(1000);
        p.setRemoveAbandonedTimeout(30000);
        p.setMinEvictableIdleTimeMillis(30000);
        p.setMinIdle(10);
        p.setLogAbandoned(true);
        p.setRemoveAbandoned(true);
        p.setJdbcInterceptors(
          "org.apache.tomcat.jdbc.pool.interceptor.ConnectionState;"+
          "org.apache.tomcat.jdbc.pool.interceptor.StatementFinalizer");
        DataSource datasource = new DataSource();
        datasource.setPoolProperties(p);
        sourceType.put(datasource, currType);
        sources.put(key, datasource);
    }


    private static Properties getPropertiesForRole(String role) throws UnknownRoleException, IOException {
        role = role.replaceAll("\\W+","_");
        if (System.getenv(role.toUpperCase() + "ROLE") != null) {
        	role = System.getenv(role.toUpperCase() + "ROLE");
        }
        String homedir = System.getenv("HOME");
        String basename = role + "_passwd";
        String fname = homedir + "/." + basename;
        File propfile = new File(fname);
        if (propfile.exists() && propfile.canRead()) {
            if (System.getenv("DEBUGPW") != null) {
                System.err.println("Opening database properties for " + role + " from " + propfile);
            }
            return readPasswd(propfile);
        }
        fname = homedir + "/" + basename;
        propfile = new File(fname);
        if (propfile.exists() && propfile.canRead()) {
            if (System.getenv("DEBUGPW") != null) {
                System.err.println("Opening database properties for " + role + " from " + propfile);
            }
            return readPasswd(propfile);
        } else {
            try {
                ClassLoader cl = ClassLoader.getSystemClassLoader();
                URL url = cl.getResource(basename);
                if (url != null) {
                    if (System.getenv("DEBUGPW") != null) {
                        System.err.println("Opening database properties for " + role + " from " + url);
                    }
                    return readPasswdStream(url.openStream());
                }
            } catch (Exception ex) {
                throw new UnknownRoleException("Couldn't find properties file for role " + role,ex);
            }
            throw new UnknownRoleException("Couldn't find properties file for role " + role);
        }         
    }

    private static Properties readPasswd(File propfile) throws IOException {
        return readPasswdStream(new FileInputStream(propfile));
    }

    private static Properties readPasswdStream(InputStream propstream) throws IOException {
        String line;
        Properties props = new Properties();
        BufferedReader reader = new BufferedReader(new InputStreamReader(propstream));        
        
        while ((line = reader.readLine()) != null) {
        	try { 
        		int p = line.indexOf('=');
        		if(p < 0) { 
        			continue;
        		}
        		String key = line.substring(0,p);
        		String value = line.substring(p+1);
        		props.setProperty(key,value);
        	} catch(RuntimeException e) { 
        		System.err.println(line);
        		throw e;
        	}
        }
        reader.close();
        return props;
    }
    /**
     * returns the type (ORACLE, MYSQL, UNKNOWN) of this
     * Connection
     */
    public static int getType(Connection c) {
        DataSource p = cxnSource.get(c);
        if (p == null) {
            return UNKNOWN;
        } else {
            return sourceType.get(p);
        }
    }
    /**
     * Returns true iff this is a connection to an Oracle DB
     */
    public static boolean isOracle(Connection c){
        return getType(c) == ORACLE;
    }
    public static boolean isMySQL(Connection c){
        return getType(c) == MYSQL;
    }
    public static boolean isPostGres(Connection c){
        return false;
    }
    public static boolean isSQLLite(Connection c){
        return false;
    }

}
