package org.seqcode.projects.akshay.chexmix.datasets;

public class BLpair {
	public BindingLocation BL1;
	public BindingLocation BL2;
	
	public BLpair(BindingLocation bl1, BindingLocation bl2) {
		this.BL1 = bl1;
		this.BL2 = bl2;
	}
	
	@Override
	public int hashCode(){
		int result =17;
		int code = (int) (this.BL1 == null ? 0 : this.BL1.hashCode());
		code += (int) (this.BL2 == null ? 0 : this.BL2.hashCode());
		result = result*37 + code;
		return result;
	}
	
	@Override
	public boolean equals(Object obj){
		if(obj == this){
			return true;
		}
		
		BLpair pair = (BLpair) obj;
		
		return this.BL1 == pair.BL1 && this.BL2 == pair.BL2;
	}
	
	

}
