package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;

public class DataCollection {

	private boolean binned=false; 
	private int columns, rows;
	private double [][][] data;
	private int binSize=-1;
	private Collection<DataSource> sources;
	private Collection<Region> regions;
	private HashMap<DataSource,Integer> source2column = new HashMap<DataSource,Integer>();
	private HashMap<Region,Integer> region2row = new HashMap<Region,Integer>();
	
	public DataCollection(Collection<DataSource> datasources, Collection<Region> regs, int binSize){
		this.binSize=binSize;
		this.binned=binSize==-1 ? false : true;
		regions=regs;
		sources=datasources;
		this.columns=datasources.size();
		this.rows=regions.size();
		if(binned)
			data = new double[rows][columns][];
		else
			data = new double[rows][columns][1];
		
		int c =0;
		for(DataSource d : sources){
			source2column.put(d, c);
			c++;
		}
		int row=0;
		for(Region r : regions){
			region2row.put(r, row);
			for(DataSource d : sources){
				if(binned)
					data[row][source2column.get(d)]=d.genDataVector(r, binSize);
				else
					data[row][source2column.get(d)][0]=d.genData(r);
			}
			row++;
		}
	}
	
	//Print to STDOUT
	public void print(){
		boolean first=true;
		for(Region r : regions){
			if(first){
				System.out.print("Region");
				for(DataSource ds : sources)
					for(int b=0; b<data[region2row.get(r)][source2column.get(ds)].length; b++)//Assuming that all rows have the same number of bins
						System.out.print("\t"+ds.getName()+"_"+b);
				System.out.print("\n");
			}
			System.out.print(r);
			for(DataSource ds : sources)
				for(int b=0; b<data[region2row.get(r)][source2column.get(ds)].length; b++)
					System.out.print("\t"+data[region2row.get(r)][source2column.get(ds)][b]);
			System.out.print("\n");
			first=false;
		}
	}
	//Print to a file
	public void print(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			boolean first=true;
			for(Region r : regions){
				if(first){
					fout.write("Region");
					for(DataSource ds : sources)
						for(int b=0; b<data[region2row.get(r)][source2column.get(ds)].length; b++)//Assuming that all rows have the same number of bins
							fout.write("\t"+ds.getName()+"_"+b);
					fout.write("\n");
				}
				fout.write(r.getLocationString());
				for(DataSource ds : sources)
					for(int b=0; b<data[region2row.get(r)][source2column.get(ds)].length; b++)
						fout.write("\t"+data[region2row.get(r)][source2column.get(ds)][b]);
				fout.write("\n");
				first=false;
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	//Remove a row	
	//This actually leaves the data structure intact, but pops the region from the list
	public void removeRow(Region r){
		regions.remove(r);
	}
	//Add rows, producing a new data structure
	public void addRows(Collection<Region> newRegs){
		int newrows = regions.size()+newRegs.size();
		double [][][] newdata;
		if(binned)
			newdata = new double[newrows][columns][];
		else
			newdata = new double[newrows][columns][1];
		
		int row=0; 
		//Old data
		for(Region r : regions){
			int oldrow = region2row.get(r);
			for(DataSource d : sources){
				newdata[row][source2column.get(d)]=data[oldrow][source2column.get(d)];
			}
			region2row.put(r, row);
			row++;
		}
		//New data
		for(Region r : newRegs){
			region2row.put(r, row);
			for(DataSource d : sources){
				if(binned)
					newdata[row][source2column.get(d)]=d.genDataVector(r, binSize);
				else
					newdata[row][source2column.get(d)][0]=d.genData(r);
			}
			row++;
		}
		
		data=newdata;
		regions.addAll(newRegs);
		rows=regions.size();
	}  
	
	//Add a single column, producing a new data structure
	public void addColumn(DataSource newSrc){
		ArrayList<DataSource> srcs = new ArrayList<DataSource>();
		srcs.add(newSrc);
		addColumns(srcs);
	}
	//Add columns, producing a new data structure
	public void addColumns(Collection<DataSource> newSrcs){
		int newcols = sources.size()+newSrcs.size();
		double [][][] newdata;
		HashMap<DataSource,Integer> newsource2column = new HashMap<DataSource,Integer>();
		if(binned)
			newdata = new double[rows][newcols][];
		else
			newdata = new double[rows][newcols][1];
		
		int col=0;
		//Old Data
		for(DataSource d : sources){
			int oldcol =source2column.get(d); 
			for(Region r : regions){
				newdata[region2row.get(r)][col]=data[region2row.get(r)][oldcol];
			}newsource2column.put(d,col);
			col++;
		}
		//New data
		for(DataSource d : newSrcs){
			newsource2column.put(d, col);
			for(Region r : regions){
				if(binned)
					newdata[region2row.get(r)][col]=d.genDataVector(r, binSize);
				else
					newdata[region2row.get(r)][col][0]=d.genData(r);
			}
			col++;
		}
		
		source2column=newsource2column;
		data=newdata;
		sources.addAll(newSrcs);
		columns = sources.size();
	}

	//Accessors
	public Collection<DataSource> getSources() {return sources;}
	public Collection<Region> getRegions() {return regions;	}
	public double[] getDataVal(Region r, DataSource s){
		if(!(region2row.containsKey(r) && source2column.containsKey(s))){
			System.err.println("Invalid query to DataCollection");
			System.exit(1);
		}
		return(data[region2row.get(r)][source2column.get(s)]);		
	}
	
	
	public void cleanup(){
		for(DataSource s : sources)
			s.cleanup();
	}
}
