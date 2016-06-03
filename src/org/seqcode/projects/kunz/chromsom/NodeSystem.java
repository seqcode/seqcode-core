package org.seqcode.projects.kunz.chromsom;
import java.util.ArrayList;

public class NodeSystem 
{
	public ArrayList<Node> nodeList;
	public ArrayList<ArrayList<Node>> nodeSystem;
	int xNode, yNode;
	
	public NodeSystem(int xNodes, int yNodes)
	{
		nodeList = new ArrayList<Node>();
		xNode = xNodes;
		yNode = yNodes;   
	    for(int i = 0; i<xNode; i++)
		{
	    	for(int m = 0; m<yNode; m++)
			{
				Node p = new Node(m,i);	
				nodeList.add(p);
			}
		}
	    nodeSystem = new ArrayList<ArrayList<Node>>();
		int nodeCount=0;
		for(int i=0; i<yNodes; i++)
		{
			nodeSystem.add(new ArrayList<Node>());
			for(int j = 0; j<xNodes; j++)
			{
				nodeSystem.get(i).add(nodeList.get(nodeCount));
				nodeCount++;
			}
		}
	}
	public int size()
	{
		return(xNode*yNode);
	}
	public Node get(int x, int y)
	{
		if(x<xNode&&y<yNode&&x>=0&&y>=0)
			return nodeSystem.get(y).get(x);
		else if(x<0 && y<0)
			return nodeSystem.get(y+yNode).get(x+xNode);
		else if(y<0)
			return nodeSystem.get(y+yNode).get(x);
		else if(x<0)
			return nodeSystem.get(y).get(x+xNode);
		else if(x>=xNode && y>=yNode)
			return nodeSystem.get(y-yNode).get(x-xNode);
		
		else if(x>=xNode && y<0)
			return nodeSystem.get(y+yNode).get(x-xNode);
		else if(x<0 && y>=yNode)
			return nodeSystem.get(y-yNode).get(x+xNode);
		
		
		else if(y>=yNode)
			return nodeSystem.get(y-yNode).get(x);
		else if(x>=xNode)
			return nodeSystem.get(y).get(x-xNode);
		
		else return null;
	}
	public Node get(int n)
	{
		return get(n/yNode, n%yNode);
	}
}
