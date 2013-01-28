package edu.psu.compbio.seqcode.gse.utils.database;

import java.util.*;

public class RoleManager {

	private Map<String,RoleConfiguration> roles;
	
	public RoleManager() { 
		roles = new HashMap<String,RoleConfiguration>();
	}
	
	public void setRoleConfiguration(String role, RoleConfiguration config) { 
		roles.put(role, config);
	}
	
}
