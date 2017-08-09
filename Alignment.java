package inversion;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Alignment  implements Comparable<Alignment> {
	public final static boolean test=Main.test;
	private String chr;
	private int refStart;
	private int refEnd;
//	int end; //skip the end information, it needs to be calculate from SA:Z
	private int readLength;
	private int mappingStart;
	private int mappingEnd;
	private String cigarValue;
	private String strand;
	private String sv;
	public  Alignment(String type,String chr,int start,String cigar,String strand){//first for length true to increase length otherwise not
		if (type.equals("xmap")){
			//?? working
		}else if(type.equals("bam")){
			this.chr=chr;
			this.refStart=start;
			this.cigarValue=cigar;
			this.strand=strand;
			this.readLength=readLength;
			this.sv="";
			this.setReadCigar(cigar);
			this.setMapping(strand, cigar);
		}
		if(test){
			System.out.println("chr="+this.chr+"\tchrStart="+this.refStart+"\tchrEnd"+this.refEnd+"\tStrand="+this.strand+"\tStart="+this.mappingStart+"\tend="+this.mappingEnd);
		}
	}
	private void setReadCigar(String cigar) {
		// TODO Auto-generated method stub
		String[] digital;
		digital=cigar.split("M|I|D|N|S|H|P|=|X");
		String[] word;
		word=cigar.split("\\d+");
//		System.out.println(Arrays.toString(digital));
//		System.out.println(Arrays.toString(word));
		int refLength=0;
		int readLength=0;
		for(int i=0;i<digital.length;i++){
			if(word[i+1].matches("M|D|N|P|=|X")){
				refLength+=Integer.parseInt(digital[i]);
			}
			if(word[i+1].matches("M|I|S|H|P|=|X")){
				readLength+=Integer.parseInt(digital[i]);
			}
			if(word[i+1].equals("D") && Integer.parseInt(digital[i])>Main.minAln ){
				int start=this.refStart+refLength-1;
				int end=this.refStart+refLength+Integer.parseInt(digital[i])-1;
				//it seems there is no read like this
				//System.out.println("Deletion="+digital[i]);
				this.sv+="DEL\t"+this.chr+"\t"+ Integer.toString(start)+"\t"+ 
						Integer.toString(end)+"\t"+ Integer.toString(end-start+1)+"\n";
			}
			if(word[i+1].equals("I") && Integer.parseInt(digital[i])>Main.minAln ){
				int start=this.refStart+refLength-1;
				int end=this.refStart+refLength;
				//System.out.println("Insertion="+digital[i]);
				this.sv+="DEL\t"+this.chr+"\t"+ Integer.toString(start)+"\t"+ 
						Integer.toString(end)+"\t"+ digital[i]+"\n";
			}
		}
		this.readLength=readLength;
		this.refEnd=this.refStart+refLength-1;
	}
	public String getChr(){
		return this.chr;
	}
	public String getCigar(){
		return this.cigarValue;
	}
	public int getRefStart(){
		return this.refStart;
	}
	public int getRefEnd(){
		return this.refEnd;
	}
	public String getStrand(){
		return this.strand;
	}
	public int getMappingStart(){
		return this.mappingStart;
	}
	public int getMappingEnd(){
		return this.mappingEnd;
	}
	public int getReadLength(){
		return this.readLength;
	}
	public String getReadSVString(){
		return this.sv;
	}
	public void setMapping(String strand, String cigar){
		String[] info;
		info=cigar.split("M|I|D|N|S|H|P|=|X");
		String start=cigar.substring(info[0].length(),info[0].length()+1);
		String end=cigar.substring(cigar.length()-1,cigar.length());
		////increase the readLength due to hard clip
		if(strand.equals("+")){
			if(start.equals("S") || start.equals("H")){
				this.mappingStart=Integer.parseInt(info[0]);
			}else{
				this.mappingStart=0;
			}
			if(end.equals("S") || end.equals("H")){
				this.mappingEnd=this.readLength-Integer.parseInt(info[info.length-1]);
			}else{
				this.mappingEnd=this.readLength;
			}
		}else if(strand.equals("-")){
			if(end.equals("S") || end.equals("H")){
				this.mappingStart=Integer.parseInt(info[info.length-1]);
			}else{
				this.mappingStart=0;
			}
			if(start.equals("S") || start.equals("H")){
				this.mappingEnd=this.readLength-Integer.parseInt(info[0]);
			}else{
				this.mappingEnd=this.readLength;
			}	
		}else{
			System.err.println("unknown strand type:"+strand);
			System.exit(1);
			return;
		}
	}
	public static void run2(String file, String args, String args2){
		File tmpFile=new File(file);
		if(!tmpFile.exists()){
			System.out.println("Input file doesn't exists. Exit.");
			System.exit(0);
		}		
		final SamReaderFactory factory =
		         SamReaderFactory.makeDefault()
	             .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
	             SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	             .validationStringency(ValidationStringency.SILENT);

		final SamReader samReader = factory.open(tmpFile);
		SAMRecordIterator iter;
		iter=samReader.iterator();
		SAMRecord sam;
		System.out.println("AlignStart"+"\t"+"AlignEnd"+"\t"+"fs"+"\t"+"ls"+"\t"+"Flags"+"\t"+"ReadName"+"\t"+"ReadLen");
		while (iter.hasNext()){
			sam=iter.next();
			int tmpInt=sam.getAlignmentEnd()-sam.getAlignmentStart();
			Cigar cigar=sam.getCigar();
			CigarElement firstElement = cigar.getCigarElement(0);
			CigarElement lastElement  = cigar.getCigarElement(cigar.numCigarElements()-1);
		}
	}
	@Override
	//compare the first start point as the key value
	public int compareTo(Alignment o) {
		return this.getMappingStart()-o.getMappingStart();
	}
}
