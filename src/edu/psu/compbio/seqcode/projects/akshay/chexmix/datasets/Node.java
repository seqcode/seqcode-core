package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.util.List;

public class Node {
	
	public boolean isleaf;
	public BindingLocation leafbl;
	public int count;
	public int[] composite;
	public Node left_child;
	public Node right_child;
	public boolean visited=false;
	
	public Node(BindingLocation leafbl) {
		this.leafbl = leafbl;
		isleaf = true;
		this.count=1;
		List<Integer> listcomp = leafbl.getConcatenatedTags(leafbl.vecpos.midpoint, leafbl.vecpos.range, "+");
		this.composite = new int[listcomp.size()];
		for(int i=0; i<listcomp.size();i++){
			this.composite[i] = listcomp.get(i);
		}
		this.left_child=null;
		this.right_child=null;
	}
	
	public Node(Node left, Node right, int[] composite) {
		this.isleaf = false;
		this.leafbl = null;
		this.count = left.count+right.count;
		this.composite = composite;
		this.left_child = left;
		this.right_child = right;
	}
	
	public static void printTree(Node root){
		while(!root.isleaf && !root.visited){
			Node left = root.left_child;
			Node right = root.right_child;
			System.out.println("left"+left.count);
			for(int i=0; i<left.composite.length; i++){
				System.out.println(i+"\t"+left.composite[i]);
			}
			System.out.println("right"+right.count);
			for(int i=0; i< right.composite.length; i++){
				System.out.println(i+"\t"+right.composite[i]);
			}
			root.visited=true;
			printTree(left);
			printTree(right);
		}
	}
	

}
