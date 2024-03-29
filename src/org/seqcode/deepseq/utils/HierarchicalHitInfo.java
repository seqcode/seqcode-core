package org.seqcode.deepseq.utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.MatrixUtils;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.genome.Genome;
import org.seqcode.math.stats.StatUtil;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;
import hdf.hdf5lib.exceptions.HDF5DataspaceInterfaceException;

/**
 * Hierarchical data structure to store the information of hits based on HDF5 project
 * The structure is:
 * 			dim1(index determined by r1 chr, strand)		dim2(index of hit info element, i.e., r1Pos, r2Chr...)			dim3(index of the hit)
 * 																		0(r1Pos)											0(100), 1(123), 2(235),...
 * 					0(/chromosome 1/-)									1(r2Chr)											0(1), 1(3), 2(2),...	
 * 																		...															...
 * 					1(/chromosome 1/+)	 	...
 * 		HitInfo		2(/chromosome 2/-)		...
 * 					3(/chromosome 2/+)		...
 * 					...
 * @author Jianyu Yang
 *
 */

public class HierarchicalHitInfo {
	
	protected Genome genome;
	
	// structure info of dataset
	protected boolean isPair;
	protected File h5;
	protected long dim1;				// dimension 1: represents the r1 chromosome and strand info
	protected long dim2;				// dimension 2: represents the each element of the hits (r1Pos, r2Chr, r2Strand, r2Pos, weight, pairMid for hitpair and start, end, weight for readhit)
	protected long dim3 = 1000000;		// dimension 3: order of each hit
	protected long[] dims;				// dim of the dataset
	protected long[] maxDims;			// extendible HDF5 dataset
	protected long[] chunkDims;			// chunk size
	protected int rank;					// number of dimensions
	protected long referenceIndex;			// which element in dimension 2 will be used to search for the hits (pairMid for hitpair and start for readhit)
	
	// HDF5 identifier
	protected long file_id = -1;
	protected long filespace_id = -1;
	protected long memspace_id = -1;
	protected long dcpl_id = -1;
	protected long dataset_id = -1;
	
	// Map used to identify the meaning of each element in dimension 2
	protected static Map<Long, String> pairID2Element = new HashMap<Long, String>();
	protected static Map<String, Long> pairElement2ID = new HashMap<String, Long>();
	protected static Map<Long, String> readID2Element = new HashMap<Long, String>();
	protected static Map<String, Long> readElement2ID = new HashMap<String, Long>();
	
	static {
		pairID2Element.put(0l, "r1Pos");
		pairID2Element.put(1l, "r2Chr");
		pairID2Element.put(2l, "r2Pos");
		pairID2Element.put(3l, "r2Strand");
		pairID2Element.put(4l, "pairWeight");
		pairID2Element.put(5l, "pairMid");
		pairID2Element.forEach((k, v) -> pairElement2ID.put(v, k));
		
		readID2Element.put(0l, "5'end");
		readID2Element.put(1l, "3'end");
		readID2Element.put(2l, "weight");
		readID2Element.forEach((k, v) -> readElement2ID.put(v, k));
	}
	
	// flag used to mark the position of newly appended hit info on each dim1
	protected int[] flag;
	// store reference array to increase the speed.
	protected double[][] ref;
	
	/**
	 * Constructor
	 * @param g			genome info to determine the dim1
	 * @param source	hdf5 file name prefix
	 * @param nElement	number of elements of the hit info which determines the dim2
	 */
	public HierarchicalHitInfo(Genome g, String source, boolean isPair) {
		this.genome = g;
		this.h5 = new File(source);
		this.isPair = isPair;
		
		// initialize the dimension information of the dataset
		dim1 = g.getChromList().size() * 2;
		dim2 = isPair? 6 : 3;
		referenceIndex = isPair ? 5 : 0;		// default reference element for pair (pairMid) or read (5'end)
		dims = new long[] {dim1, dim2, dim3};
		maxDims = new long[] {dim1, dim2, HDF5Constants.H5S_UNLIMITED};
		chunkDims = new long[] {1, dim2, 100};
		rank = 3;
		
		flag = new int[(int)dim1];
		ref = new double[(int)dim1][];
		
	}
	
	// Setters
	public void setReference(long referenceID) { referenceIndex = referenceID; }
	public void setReference(String referenceString) { referenceIndex = isPair ? pairElement2ID.get(referenceString) : readElement2ID.get(referenceString); }
	
	// Accessors
	public int getLength(String chr, int strand) { return flag[convertIndex(chr, strand)];}
	public static Set<Long> getReadElementIDList() {return readID2Element.keySet();}
	public static Set<Long> getPairElementIDList() {return pairID2Element.keySet();}
	public long getReferenceIndex() {return referenceIndex;}
	public double getReference(String chr, int strand, int order) throws Exception{
		return ref[convertIndex(chr, strand)][order];
	}
	
	public double[] getReference(String chr, int strand, int startOrder, int endOrder) throws Exception {
		return Arrays.copyOfRange(ref[convertIndex(chr, strand)], startOrder, endOrder);
	}

	/**
	 * get specific element in the dataset
	 * @param chr
	 * @param strand
	 * @param elementID
	 * @param order
	 * @return
	 */
	public double getElement(String chr, int strand, long elementID, int order) throws Exception {
		int index = convertIndex(chr, strand);
		if(flag[index]<order) {throw new IndexOutOfBoundsException();}
		long[] start = { index, elementID, order };
		long[] count = { 1, 1, 1};
		return getElement(start, count)[0][0];
	}
	
	/**
	 * get specific element in the dataset
	 * @param chr
	 * @param strand
	 * @param elementString
	 * @param order
	 * @return
	 */
	public double getElement(String chr, int strand, String elementString, int order) throws Exception {
		return getElement(chr, strand, isPair ? pairElement2ID.get(elementString) : readElement2ID.get(elementString), order);
	}
	
	/**
	 * get indice of the element in the dataset
	 * @param chr
	 * @param strand
	 * @param elementID
	 * @param startOrder
	 * @param endOrder
	 * @return
	 */
	public double[] getElement(String chr, int strand, long elementID, int startOrder, int endOrder) throws Exception {
		if(endOrder<=startOrder)
			return null;
		int index = convertIndex(chr, strand);
		if(flag[index]<endOrder) {throw new IndexOutOfBoundsException();}
		long[] start = { index, elementID, startOrder };
		long[] count = { 1, 1, endOrder-startOrder};
		return getElement(start, count)[0];
	}
	
	public double[] getElement(String chr, int strand, String elementString, int startOrder, int endOrder) throws Exception {
		return getElement(chr, strand, isPair ? pairElement2ID.get(elementString) : readElement2ID.get(elementString), startOrder, endOrder);
	}
	
	public double[] getElement(String chr, int strand, String elementString) {
		return getElement(chr, strand, isPair ? pairElement2ID.get(elementString) : readElement2ID.get(elementString));
	}
	
	public double[] getElement(String chr, int strand, long elementID) {
		int index = convertIndex(chr, strand);
		long[] start = { index, elementID, 0};
		long[] count = { 1, 1, flag[index]};
		return getElement(start, count)[0];
	}
	
	public double[][] getElement(String chr, int strand, int startOrder, int endOrder) throws Exception {
		if(endOrder<=startOrder)
			return null;
		int index = convertIndex(chr, strand);
		if(flag[index]<endOrder) {throw new IndexOutOfBoundsException();}
		long[] start = { index, 0, startOrder};
		long[] count = { 1, isPair ? pairElement2ID.size() : readElement2ID.size(), endOrder-startOrder};
		return getElement(start, count);
	}
	/**
	 * Get a region of data in the dataset
	 * I think we won't need to report indice of data across several chromosome/strand dimension, so only two dimensions will be reported
	 * @param start
	 * @param count
	 * @return
	 */
	public double[][] getElement(long[] start, long[] count) {
		double[][] elements = new double[(int)count[1]][(int)count[2]];
		try {
			// get filespace_id
			filespace_id = H5.H5Dget_space(dataset_id);
			
			if (filespace_id >= 0) {
				H5.H5Sselect_hyperslab(filespace_id, HDF5Constants.H5S_SELECT_SET, start, null, count, null);
				
				// read the corresponding hit reference
				memspace_id = H5.H5Screate_simple(rank, count, null);
				if (dataset_id >= 0)
					H5.H5Dread(dataset_id, HDF5Constants.H5T_NATIVE_DOUBLE, memspace_id, 
							filespace_id, HDF5Constants.H5P_DEFAULT, elements);
				
				// End access
				try {
					H5.H5Sclose(filespace_id);
					H5.H5Sclose(memspace_id);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		return elements;
	}
	
	public StrandedPair getPair(String chr, int strand, int order) throws Exception {
		return getPair(chr, strand, order, order+1).get(0);
	}
	
	public List<StrandedPair> getPair(String chr, int strand, int startOrder, int endOrder) throws Exception {
		if(!isPair)
			throw new Exception("Extract hitpair info from a readHit dataset");
		if(endOrder<=startOrder)
			return null;
		double[][] hitInfo = getElement(chr, strand, startOrder, endOrder);
		
		// transpose the hitInfo
		long start = System.currentTimeMillis();
		double[][] transposedHitInfo = transposeMatrix(hitInfo);
		
		// convert each row of the hitInfo to strandedPair then return
		List<StrandedPair> pairs = new ArrayList<StrandedPair>();
		for(int i=0; i<transposedHitInfo.length; i++)
			pairs.add(convertPair(chr, strand, transposedHitInfo[i][0], transposedHitInfo[i][1], transposedHitInfo[i][2],
					transposedHitInfo[i][3], transposedHitInfo[i][4]));
		
		// Free memory
		// TODO: I don't think it is a right way to free memory, but I still haven't found a better way to make sure
		// the hitInfo and transposedHitInfo will be recycled instantly
		hitInfo = null;
		transposedHitInfo = null;
		
		return pairs;
	}
	
	public StrandedBaseCount getBase(String chr, int strand, int order) throws Exception {
		return getBase(chr, strand, order, order+1).get(0);
	}
	
	public List<StrandedBaseCount> getBase(String chr, int strand, int startOrder, int endOrder) throws Exception{
		if(isPair) {
			throw new Exception("Extract readhit info from a hitpair dataset");
		}
		if(endOrder<=startOrder)
			return null;
		double[][] hitInfo= getElement(chr, strand, startOrder, endOrder);
		
		// transpose the hitInfo
		double[][] transposedHitInfo = transposeMatrix(hitInfo);
		
		// convert each row of the hitInfo to strandedBaseCount then return
		List<StrandedBaseCount> bases = new ArrayList<StrandedBaseCount>();
		for(int i=0; i<transposedHitInfo.length; i++)
			bases.add(convertBase(chr, strand, transposedHitInfo[i][0], transposedHitInfo[i][1], transposedHitInfo[i][2]));
		
		// Free memory
		// TODO: I don't think it is the right way to free memory, but I still haven't found a better way to make sure
		// the hitInfo and transposedHitInfo will be recycled instantly
		hitInfo = null;
		transposedHitInfo = null;
		
		return bases;
	}
	
	public void updateElement(String chr, int strand, String elementString, double[] element) throws Exception{
		updateElement(chr, strand, isPair ? pairElement2ID.get(elementString) : readElement2ID.get(elementString), element);
	}
	
	public void updateElement(String chr, int strand, Long elementID, double[] element) throws Exception{
		if(element.length != getLength(chr, strand))
			throw new Exception("Length of input element is unequal to the flag Expected " + getLength(chr, strand) + "\t Get: " + element.length);
		
		long[] start = {convertIndex(chr, strand), elementID, 0};
		long[] count = {1, 1, flag[convertIndex(chr, strand)]};
		try {
			// get filespace_id
			filespace_id = H5.H5Dget_space(dataset_id);
		
			H5.H5Sselect_hyperslab(filespace_id, HDF5Constants.H5S_SELECT_SET, start, null, count, null);
			
			// write the element to the dataset
			memspace_id = H5.H5Screate_simple(rank, new long[] {1, 1, flag[convertIndex(chr, strand)]}, null);
			if (dataset_id >= 0)
				H5.H5Dwrite(dataset_id, HDF5Constants.H5T_NATIVE_DOUBLE, memspace_id, 
						filespace_id, HDF5Constants.H5P_DEFAULT, element);
			
			// End access to the memory space to release the resource
			H5.H5Sclose(filespace_id);
			H5.H5Sclose(memspace_id);
					
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Sort all elements according to the reference element
	 */
	public void sortByReference() {
		for (String chrom: genome.getChromList())
			for (int strand=0; strand<2; strand++) { if(getLength(chrom, strand) > 0) {
				// 1. load the reference element
				double[] reference = getElement(chrom, strand, referenceIndex);
				
				// 2. sort all elements according to the reference then write back to the dataset
				int[] inds = StatUtil.findSort(reference);
				ref[convertIndex(chrom, strand)] = reference;
				for(long elementID: isPair? pairID2Element.keySet() : readID2Element.keySet()) {
					if(elementID != referenceIndex) {
						double[] sortElement = permute(getElement(chrom, strand, elementID), inds);
						try {
							updateElement(chrom, strand, elementID, sortElement);
						} catch (Exception e) {
							e.printStackTrace();
							System.exit(1);
						}
						sortElement = null;
					}
				}
				
				// 3. write the reference element back to the dataset finally
				try {
					updateElement(chrom, strand, referenceIndex, reference);
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(1);
				}
				
				reference = null;
				inds = null;
			}
		}
	}
	
	protected StrandedPair convertPair(String chr, int strand, double r1Pos, double r2Chr, double r2Pos, double r2Strand, double pairWeight) {
		return new StrandedPair(genome, genome.getChromID(chr), (int)r1Pos, strand==0? '+' : '-',
				(int)r2Chr, (int)r2Pos, r2Strand==0 ? '+' : '-', (int)pairWeight);
	}
	
	protected StrandedBaseCount convertBase(String chr, int strand, double end5, double end3, double weight) {
		return new StrandedBaseCount(strand==0 ? '+' : '-', (int)(strand==0? end5 : end3), (float)weight);
	}
	
	/**
	 * Initialize the structure of HDF5 file
	 */
	public void initializeHDF5() {
		// Create a new file using default properties
		try {
			file_id =  H5.H5Fcreate(h5.getAbsolutePath(), HDF5Constants.H5F_ACC_TRUNC, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Create dataspace. Setting maximum size as unlimited to enable extendible dataset
		try {
			filespace_id = H5.H5Screate_simple(rank, dims, maxDims);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Create the dataset creation property list
		try {
			dcpl_id = H5.H5Pcreate(HDF5Constants.H5P_DATASET_CREATE);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Set the chunk size
		try {
			if (dcpl_id >= 0)
				H5.H5Pset_chunk(dcpl_id, rank, chunkDims);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		
		// Create the chunked dataset
		try {
			if ((file_id >= 0) && (filespace_id >= 0) && (dcpl_id >= 0 ))
				dataset_id = H5.H5Dcreate(file_id, "HitInfo", HDF5Constants.H5T_NATIVE_DOUBLE, filespace_id, 
						HDF5Constants.H5P_DEFAULT, dcpl_id, HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		
		// End access to the resource
        try {
            if (dcpl_id >= 0)
                H5.H5Pclose(dcpl_id);
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        try {
            if (filespace_id >= 0)
                H5.H5Sclose(filespace_id);
        }
        catch (Exception e) {
            e.printStackTrace();
        }

	}
	
	/**
	 * Open an existing HDF5HitCache info
	 */
	public void openDataSet() {
		// Open an existing file using default properties
		try {
			file_id =  H5.H5Fopen(h5.getAbsolutePath(), HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Get dataset id
		try {
			if(file_id>=0)
				dataset_id = H5.H5Dopen(file_id, "HitInfo", HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Iterate each chromosome/strand to get the length of the info (Use 5'end for reads or r1Pos for hitpair because these values are always > 0)
		try {
			filespace_id = H5.H5Dget_space(dataset_id);
			H5.H5Sget_simple_extent_dims(filespace_id, dims, maxDims);
			H5.H5Sclose(filespace_id);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		for(String chr: genome.getChromList())
			for(int strand=0; strand<2; strand++) {
				flag[convertIndex(chr, strand)] = Integer.MAX_VALUE;
				int endIndex = 0;
				try {
					if(getElement(chr, strand, isPair ? "r1Pos" : "5'end", (int)(dims[2]-1)) > 0) {
						flag[convertIndex(chr, strand)] = (int)dims[2];
					} else {
						// First find a roughly region of the end index
						while(true) {
							if(getElement(chr, strand, isPair ? "r1Pos" : "5'end", endIndex) > 0) {
								endIndex+=10000;
								if(endIndex >= dims[2])
									endIndex = (int)(dims[2] - 1);
							} else
								break;
						}
						// Then find the accurate end index
						while(true) {
							if(getElement(chr, strand, isPair ? "r1Pos" : "5'end", endIndex) == 0) {
								endIndex--;
								if(endIndex<0)
									break;
							}
							else
								break;
						}
						flag[convertIndex(chr, strand)] = endIndex+1;
					}
					System.err.println("Detected number of hits on chromosome " + chr + " strand " + strand + " : " + getLength(chr, strand));
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
	}
	
	/**
	 * close the dataset
	 */
	public void closeDataset() {
        // End access to the dataset and release resources used by it.
        try {
            if (dataset_id >= 0)
                H5.H5Dclose(dataset_id);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
	}
	
	/**
	 * close the file
	 */
	public void closeFile() {
        try {
            if (file_id >= 0)
                H5.H5Fclose(file_id);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
	}
	
	public void deleteFile() {
		if(h5.delete())
			System.err.println("Successfully delete HDF5 file: " + h5.getAbsolutePath());
	}
	
	/**
	 * extend the dataset by 1000000 columns
	 */
	public void extendDataset() {
		extendDataset(10000);
	}
	
	public void extendDataset(long size) {
		try {		
			// extend the dataset
			dim3 += size;
			if (dataset_id >= 0)
				H5.H5Dset_extent(dataset_id, new long[] {dim1, dim2, dim3});
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public int convertIndex(String chr, int strand) {
		return genome.getChromID(chr)*2 + strand;
	}
	
	public void appendHit(String chr, int strand, double[] hit) throws Exception {
		// check the size of hit
		if(hit.length != dim2)
			throw new Exception("Unequal dimension of hit info! Excepted: " + dim2 + " Get: " + hit.length);
		// select the part written by the hit info and write the info
		int index = convertIndex(chr, strand);
		long[] start = { index, 0, flag[index]++ };
		long[] count = { 1, dim2, 1};
		
		// extend the dataset if flag exceed the size of dimension 3
		while (flag[index] > dim3)
			extendDataset();
				
		try {
			// get filespace_id
			filespace_id = H5.H5Dget_space(dataset_id);
			
			 
			// Select the region used to append hit
			if (filespace_id >= 0)
				H5.H5Sselect_hyperslab(filespace_id, HDF5Constants.H5S_SELECT_SET, start, null, count, null);
				
			// Create the memory space to append the hit
			memspace_id = H5.H5Screate_simple(rank, new long[] {1, dim2, 1}, null);
			
			// write the data into the dataset
			if (dataset_id >= 0 && filespace_id >= 0 && memspace_id >= 0)
				H5.H5Dwrite(dataset_id, HDF5Constants.H5T_NATIVE_DOUBLE, memspace_id, 
						filespace_id, HDF5Constants.H5P_DEFAULT, hit);
			
			// End access to the memory space to release the resource
			H5.H5Sclose(filespace_id);
			H5.H5Sclose(memspace_id);
			
		}  catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void appendArray(String chr, int strand, double[][] hits) throws Exception {
		// check the size of the hit
		if(hits.length != dim2)
			throw new Exception("Unequal dimension of hit info! Excepted: " + dim2 + " Get: " + hits.length);
		// select the part written by the hit info and write the info
		int index = convertIndex(chr, strand);
		long[] start = { index, 0, flag[index] };
		long[] count = { 1, dim2, hits[0].length};
		flag[index] += hits[0].length;
		
		// extend the dataset if flag exceed the size of dimension 3
		while (flag[index] > dim3)
			extendDataset();
		
		// append the array into the dataset
		try {
			// get filespace_id
			filespace_id = H5.H5Dget_space(dataset_id);
			
			// select the region to append the array
			H5.H5Sselect_hyperslab(filespace_id, HDF5Constants.H5S_SELECT_SET, start, null, count, null);
			
			// Create the memory space to append the hit
			memspace_id = H5.H5Screate_simple(rank, new long[] {1, dim2, hits[0].length}, null);
			
			// write the data into the dataset
			if (dataset_id >= 0 && filespace_id >= 0 && memspace_id >= 0)
				H5.H5Dwrite(dataset_id, HDF5Constants.H5T_NATIVE_DOUBLE, memspace_id, 
						filespace_id, HDF5Constants.H5P_DEFAULT, hits);
			
			// End access to the memory space to release the resource
			H5.H5Sclose(filespace_id);
			H5.H5Sclose(memspace_id);
			
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
	}
	
	public double[] permute(double[] array, int[] inds) {
		double[] returnArray = new double[array.length];
		int index = 0;
		for(int i: inds)
			returnArray[index++] = array[i];
		return returnArray;
	}
	
	public static double[][] transposeMatrix(double[][] matrix){
	    int m = matrix.length;
	    int n = matrix[0].length;

	    double[][] transposedMatrix = new double[n][m];

	    for(int x = 0; x < n; x++) {
	        for(int y = 0; y < m; y++) {
	            transposedMatrix[x][y] = matrix[y][x];
	        }
	    }

	    return transposedMatrix;
	}
}




























