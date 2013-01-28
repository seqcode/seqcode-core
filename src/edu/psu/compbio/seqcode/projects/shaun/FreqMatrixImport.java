package edu.psu.compbio.seqcode.projects.shaun;
import java.io.*;
import java.util.*;
import java.sql.*;
import java.text.ParseException;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrixImport;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;
import edu.psu.compbio.seqcode.gse.utils.io.motifs.BackgroundModelIO;

/* Loads frequency and count matrix motifs to the db. 
 * The motifs are first converted to log-odds weightmatrices using a background model.
 * The background model file can be generated using MakeBackgroundModel or left unspecified, in which case the mononucleotide background is taken from the genome.   

Usage is similar to WeightMatrixImport:
java edu.psu.compbio.seqcode.gse.shaun.FreqMatrixImport --species 'Saccharomyces cerevisiae' --fmname HSF1 --fmversion TRANSFAC11.3 --fmtype TRANFSAC --fmfile HSF1.motif --backfile yeast.back

*/

public class FreqMatrixImport {

	private static WeightMatrixImport wmimp = new WeightMatrixImport();
	private static MarkovBackgroundModel back = null;
	private static int[] indices = { 'A', 'C', 'G', 'T' };
	private static int MAX_MOTIF_LEN = 200;
	private static float SCALE_FACTOR = (float) 0.1;
	
    public static void main(String args[]) throws IOException, ParseException, NotFoundException {
        String species = null;
        String fmname = null, fmversion = null, fmtype = null;
        String fmfile = null; String backfile = null; String genome=null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) { 
                species = args[++i];
                if (species.indexOf(';') != -1) {
                    String[] pieces = species.split(";");
                    species = pieces[0];
                }
            }
            if (args[i].equals("--fmname")) {
                fmname = args[++i];
                if (fmname.indexOf(';') != -1) {
                    String[] pieces = fmname.split(";");
                    fmname = pieces[0];
                    fmversion = pieces[1];
                    if (pieces.length >= 3) {
                        fmtype = pieces[2];
                    }
                }
            }
            if (args[i].equals("--fmversion")) {
                fmversion = args[++i];
            }
            if (args[i].equals("--fmtype")) {
                fmtype = args[++i];
            }
            if (args[i].equals("--backfile")) {
                backfile = args[++i];
            }if (args[i].equals("--genome")) {
                genome = args[++i];
            }
            if (args[i].equals("--") ||
                args[i].equals("--fmfile")) {
                fmfile = args[++i];
            }
        }

        if (species == null) {
            System.err.println("Must supply a --species"); System.exit(1);
        }
        if (fmfile == null) {
            System.err.println("Must supply a --wmfile"); System.exit(1);
        }
        if(backfile==null && genome!=null){
          back = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(Organism.findGenome(genome)));        	
        }else if(backfile==null){
        	System.err.println("Must supply --backfile or --genome"); System.exit(1);
        }else{
        	back = BackgroundModelIO.parseMarkovBackgroundModel(backfile, Organism.findGenome(genome));        	
        }
        try {
            if(fmname==null) { 
                insertMultiFMFromFile(species,fmversion, fmtype,fmfile);
            } else { 
                if (fmversion == null) {
                    System.err.println("Must supply a --wmversion"); System.exit(1);
                }
                insertFMFromFile(species,fmname,fmversion,fmtype,fmfile);
            }
        } catch (SQLException ex) {
            ex.printStackTrace();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
            System.err.println("Must supply a valid species and genome");
        } catch (UnknownRoleException ex ){
            ex.printStackTrace();
            System.err.println("Couldn't connect to role annotations");
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
            System.err.println("Couldn't find the input file");
        } catch (ParseException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    public void setBackground(MarkovBackgroundModel b){back=b;}
    
    //Unlike the WeigthMatrixImport version, this one actually loads more than one matrix!
    public static LinkedList<WeightMatrix> readTransfacMatrices(String wmfile, String version) throws IOException {
        LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
        BufferedReader br = new BufferedReader(new FileReader(new File(wmfile)));
        String line;        
        WeightMatrix matrix = null;
        int motifCount=0;
        Vector<float[][]> arrays = new Vector<float[][]>();
        Vector<Integer> arrayLens = new Vector<Integer>();
        Vector<String> names = new Vector<String>();
        Vector<String> versions=new Vector<String>();

        //Read in Transfac format first
        boolean nameLoaded=false;
        int matLen=0;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) {
            	String[] pieces = line.split("\\s+");
                if(pieces[0].equals("DE")){
                	names.add(pieces[1]);
                	if(pieces.length>=3){
                		String v_string =pieces[2];
                		if(pieces.length>=4){for(int v=3; v<pieces.length; v++){v_string =v_string+","+pieces[v];}}
                		
                		if(version!=null){versions.add(new String(v_string+","+version));}
                		else{versions.add(v_string);}
                	}else{
                		versions.add(version);
                	}
                    nameLoaded=true;
                    arrays.add(new float[MAX_MOTIF_LEN][4]);
                    matLen=0;
                }else if(pieces[0].equals("XX")){
                	arrayLens.add(matLen);
                	motifCount++;
                }else if(nameLoaded && (pieces.length==5 || pieces.length==6)){ 
                	//Load the matrix
                	for(int i = 1; i <=4 ; i++) { 
                        arrays.get(motifCount)[matLen][i-1] = Float.parseFloat(pieces[i]);
                    }
                    matLen++;
                }
            }
        }
        for(int m = 0; m<motifCount; m++){
        	//Make a new WeightMatrix
            matrix = new WeightMatrix(arrayLens.get(m));
            matrix.name=names.get(m);
            matrix.version=versions.get(m);

            //Convert the freq matrix to a weight matrix
            for(int i = 0; i < arrayLens.get(m); i++) { 
            	float ttl=0;
            	for(int j=0; j<4; j++){ttl += arrays.get(m)[i][j];}
            	float currScale = SCALE_FACTOR*ttl; 
        		for(int j = 0; j < 4; j++) {
        			matrix.matrix[i][indices[j]] = (float)(Math.log((((arrays.get(m)[i][j] + (currScale*back.getMarkovProb(j, 1)))/(ttl+currScale))/back.getMarkovProb(j, 1)))/Math.log(2));
                }
            }
            matrices.add(matrix);
           System.err.println("Added \"" + matrix.name + "\"\t"+matrix.version);
	   System.err.println("Max score = "+matrix.getMaxScore());
           System.err.println(WeightMatrix.printMatrix(matrix));
        }

        br.close();
        return matrices;
    }
    public static LinkedList<WeightMatrix> readTransfacMatrices(String wmfile, String version, float pseudocountTotal){
    	SCALE_FACTOR = pseudocountTotal;
    	LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
    	try {
			matrices = readTransfacMatrices(wmfile, version);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}return(matrices);
    }
    public static LinkedList<WeightMatrix> readTransfacMatrices(String wmfile){
    	LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
    	try {
			matrices = readTransfacMatrices(wmfile, "");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}return(matrices);
    }
    
    public static Set<Integer> insertMultiFMFromFile(String species, String version,String fmtype, String fmfile) 
        throws IOException, SQLException, NotFoundException { 
    	LinkedList<WeightMatrix> matrices=null;
        HashSet<Integer> ids = new HashSet<Integer>();
        
        if(fmtype.matches(".*TRANSFAC.*")) { 
            matrices = readTransfacMatrices(fmfile, version);
        }else if (fmtype.matches(".*SEQ.*")) { 
            try {
				matrices = readAlignedSequenceMatrices(fmfile);
			} catch (ParseException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        } else {
            System.err.println("Didn't see a program I recognize in the type.  defaulting to TRANSFAC");
            matrices = readTransfacMatrices(fmfile);
        }

        if(matrices!=null){
            for(WeightMatrix matrix : matrices) { 
                matrix.type = fmtype;
                matrix.species = species;
                matrix.setLogOdds();
                ids.add(wmimp.insertMatrixIntoDB(matrix));                
            }
        }System.out.println("Matrices Loaded: "+matrices.size());
        return ids;
    }

    /* reads a freq matrix from the specified file,
       inserts it into the database with the specified name and
       version, and returns its dbid.  This means that the file may contain
       only a single weight matrix*/
    public static int insertFMFromFile(String species,
                                       String fmname,
                                       String fmversion,
                                       String fmtype,
                                       String fmfile) throws SQLException, NotFoundException, UnknownRoleException, FileNotFoundException, ParseException, IOException {
    	LinkedList<WeightMatrix> matrices;
    	
        if (fmtype.matches(".*TRANSFAC.*")) {
            matrices = readTransfacMatrices(fmfile);
        } else if (fmtype.matches(".*SEQ.*")) { 
            matrices = readAlignedSequenceMatrices(fmfile);
        } else {
            System.err.println("Didn't see a program I recognize in the type.  defaulting to TRANSFAC");
            matrices = readTransfacMatrices(fmfile);
        }
        matrices.get(0).name = fmname;
        matrices.get(0).version = fmversion;
        matrices.get(0).type = fmtype;
        matrices.get(0).species = species;
        matrices.get(0).setLogOdds();
        return wmimp.insertMatrixIntoDB(matrices.get(0));
    }
    
    /**
     * Constructs a matrix from a set of strings.  The strings must all have the same length.
     * The WeightMatrix returned has the frequencies of the bases at each position given.
     * 
     * @param strings
     * @return
     * @throws IOException
     * @throws ParseException
     */
    public static WeightMatrix buildAlignedSequenceMatrix(Collection<String> strings) throws ParseException {
        WeightMatrix wm = null;
        
        int[] counts = null;
        for(String line : strings) { 
            line = line.trim().toUpperCase();
            if(line.length() > 0) { 
                if(wm == null) { 
                    wm = new WeightMatrix(line.length());
                    counts = new int[line.length()];
                    for(int i = 0; i < wm.length(); i++) {
                        counts[i] = 0;
                        for(int j = 0; j < wm.matrix[i].length; j++) { 
                            wm.matrix[i][j] = (float)0.0;
                        }
                    }
                }
                
                if(line.length() != wm.length()) { 
                    throw new ParseException("Line \"" + line + "\" was of uneven length (" + 
                            wm.length() + ")", 0);
                }
                
                for(int i = 0; i < line.length(); i++) {
                    char c = line.charAt(i);
                    if(c != 'N' && c != '-') { 
                        wm.matrix[i][c] += (float)1.0;
                        counts[i] += 1;
                    }
                }
            } 
        }
        
        for(int i = 0; wm != null && i < wm.length(); i++) { 
            if(counts[i] > 0) { 
                for(int j = 0; j < wm.matrix[i].length; j++) {
                    wm.matrix[i][j] /= (float)counts[i];
                }
            } else { 
                for(int j = 0; j < wm.matrix[i].length; j++) { 
                    wm.matrix[i][j] = (float)1.0 / (float)(wm.matrix[i].length);
                }
            }
        }
        
        return wm;
    }

    /**
     * Reads a matrix from a text-file.  The file is assumed to have one sequence per line, 
     * and all lines must be the same length.  The WeightMatrix returned has the frequencies 
     * of the bases at each position given.
     * 
     * @param fname
     * @return
     * @throws IOException
     * @throws ParseException
     */
    public static LinkedList<WeightMatrix> readAlignedSequenceMatrices(String fname) throws IOException,ParseException {
    	LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
    	int motifCount=0;
    	
    	LinkedList<String> strings = new LinkedList<String>();
        String line, currName="Motif";
        BufferedReader br = new BufferedReader(new FileReader(new File(fname)));
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.contains(">")){
            	line.replaceAll(">", "");
            	currName = line;
            	if(strings.size()>0){
            		matrices.add(buildAlignedSequenceMatrix(strings));
            		matrices.get(motifCount).name = currName;
            	}
            	strings = new LinkedList<String>();
            	motifCount++;
            }else if(line.length() > 0) { 
                strings.addLast(line);
            }
        }
        br.close();

        if(strings.size()>0){
    		matrices.add(buildAlignedSequenceMatrix(strings));
    		matrices.get(motifCount).name = currName;
    	}
        
        return(matrices);
    }

}

