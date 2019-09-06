package org.seqcode.deepseq.utils;

import java.util.HashMap;
import java.util.Map;
import java.util.function.LongPredicate;

import org.seqcode.deepseq.HitPair;
import org.seqcode.genome.Genome;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;
import htsjdk.tribble.index.Block;

/**
 * Hierarchical data structure to store the information of hits based on HDF5 project
 * The structure is:
 * 											[ 0, 1, 1, 3]
 * 					0(/chromosome 1/-)		[ 3, 5, 6, 2]
 * 											[ 5, 2, 4, 6]
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
	protected String filename;
	protected long dim1;
	protected long dim2;
	protected long dim3 = 10000;
	protected long[] dims;
	protected long[] maxDims;
	protected long[] chunkDims;
	protected int rank;
	
	// HDF5 identifier
	protected long file_id = -1;
	protected long filespace_id = -1;
	protected long memspace_id = -1;
	protected long dcpl_id = -1;
	protected long dataspace_id = -1;
	protected long dataset_id = -1;
	
	// flag used to mark the position of newly appended hit info on each dim1
	protected int[] flag;
	
	public HierarchicalHitInfo(Genome g, String source, int nElement) {
		this.genome = g;
		this.filename = source + ".h5";
		
		dim1 = g.getChromList().size() * 2;
		dim2 = nElement;
		dims = new long[] {dim1, dim2, dim3};
		maxDims = new long[] {dim1, dim2, HDF5Constants.H5S_UNLIMITED};
		chunkDims = new long[] {1, dim2, 10000};
		rank = 3;
		
		flag = new int[(int)dim1];
	}
	
	/**
	 * Initialize the structure of HDF5 file
	 */
	public void initializeHDF5() {
		// Create a new file using default properties
		try {
			file_id =  H5.H5Fcreate(filename, HDF5Constants.H5F_ACC_TRUNC, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
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
				dataset_id = H5.H5Dcreate(file_id, "HitInfo", HDF5Constants.H5T_IEEE_F32BE, filespace_id, 
						HDF5Constants.H5P_DEFAULT, dcpl_id, HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		
		// Create the memory space to append the hit
		try {
			memspace_id = H5.H5Screate_simple(rank, new long[] {1, dim2, 1}, null);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}

	}
	
	/**
	 * close the dataset
	 */
	public void closeDataset() {
        // End access to the dataset and release resources used by it.
        try {
            if (dcpl_id >= 0)
                H5.H5Pclose(dcpl_id);
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        try {
            if (dataset_id >= 0)
                H5.H5Dclose(dataset_id);
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
		try {
			if (filespace_id >= 0) {
				H5.H5Sselect_hyperslab(filespace_id, HDF5Constants.H5S_SELECT_SET, start, null, count, null);
				
				// write the data into the dataset
				if (dataset_id >= 0)
					H5.H5Dwrite(dataset_id, HDF5Constants.H5T_IEEE_F32BE, memspace_id, 
							filespace_id, HDF5Constants.H5P_DEFAULT, hit);
			}
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
	}
	
}




























