package org.seqcode.projects.kunz.chromsom;

import java.util.ArrayList;

public class MiniSystem 
{
	public ArrayList<MiniNode> miniList;
	public ArrayList<ArrayList<MiniNode>> miniSystem;
	int xNode, yNode;
	
	public MiniSystem(ArrayList<MiniNode> nds,int xNodes, int yNodes)
	{
		miniList = nds;
		xNode = xNodes;
		yNode = yNodes;   
	    miniSystem = new ArrayList<ArrayList<MiniNode>>();
		int nodeCount=0;
		for(int i=0; i<yNodes; i++)
		{
			miniSystem.add(new ArrayList<MiniNode>());
			for(int j = 0; j<xNodes; j++)
			{
				miniSystem.get(i).add(miniList.get(nodeCount));
				nodeCount++;
			}
		}
	}
	public int size()
	{
		return(xNode*yNode);
	}
	public MiniNode get(int x, int y)
	{
		if(x<xNode&&y<yNode&&x>=0&&y>=0)
			return miniSystem.get(y).get(x);
		else if(x<0 && y<0)
			return miniSystem.get(y+yNode).get(x+xNode);
		else if(y<0)
			return miniSystem.get(y+yNode).get(x);
		else if(x<0)
			return miniSystem.get(y).get(x+xNode);
		else if(x>=xNode && y>=yNode)
			return miniSystem.get(y-yNode).get(x-xNode);
		else if(y>=yNode)
			return miniSystem.get(y-yNode).get(x);
		else if(x>=xNode)
			return miniSystem.get(y).get(x-xNode);
		
		else return null;
	}
	public MiniNode get(int n)
	{
		return get(n/yNode, n%yNode);
	}
}
