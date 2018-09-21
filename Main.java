package inversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
//for java 1.7
//import java.time.LocalDateTime;
//import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

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



// focus on inversion
// genotype the read base on aligned and non-align(remove the non-information read)
// check on the validated data
// min border max border
public class Main {
	//global setting for the program
	public final static boolean test=false;//STDOUT the testing information if true
	public static int minAln;
	public static int min;
	public static int max;
	public static int gq;
	public final static int clipSize=5;
	public static int mergedClipRead=10;
	public static double variantRate=0.01;
	public static double version=1.21;
	public static final int mergedSize=2000;//sv merged size
	public static int minimumClipRead=3;//short clip read record>=3
	//end of global setting
	//for all candidate long sv signal
	//no more global data!
	//public static ArrayList<SV> sv=new ArrayList<SV>();//global sv set
	// for all long read start end
	//public static HashMap<String, startEnd[]> depth=new HashMap<String,startEnd[]>();
	//for all short read start end
	//public static HashMap<String, startEnd[]> depthShort=new HashMap<String, startEnd[]>();
	//for all clip read information
	//public static ClipRead[] clip;//they are sorted in the same chromosome
	//entrance point
	public static void main(String[] args) throws IOException{
		getArgs tmp=new getArgs(args);
		tmp.addUsage("Program function: Read a SE bam file and get the inversion\n");
		tmp.addUsage("Version:\t"+version+"\n");
		tmp.addUsage("--output[String] file to write\n");
		tmp.addUsage("--input[String] file to read\n");
		tmp.addUsage("optional:\n");
		//tmp.addUsage("--inputType[String] bam or xmap(BioNano). Default[bam]\n");
		tmp.addUsage("--region[String] Specify the region for running.\n                 Such as chr9:1-1000 OR chr9 OR all Default[all]\n");
		tmp.addUsage("--minAln[int] minimum size for Alignment & Inv. Default[500]\n");
		tmp.addUsage("--IRdatabase[String] An inverted repeat file for the reference in bed format. Default[none]\n");
		//tmp.addUsage("--shortBam[String] An short bam file for the same sample to the same reference. Default[none]\n");
		tmp.addUsage("--min[int] minimum size of inversion. Default[500]\n");
		tmp.addUsage("--max[int] maximum size of inversion. Default[10000]\n");
		tmp.addUsage("For example: java -jar npInv.jar --input sample.bam --output sample.VCF\n");
		String output=tmp.getString("output", true, "output");
		String inputType=tmp.getString("inputType", false, "bam");
		String input=tmp.getString("input", true, "xx");
		String region=tmp.getString("region", false, "all");
		String ir=tmp.getString("IRdatabase", false, "none");
		String shortBam=tmp.getString("shortBam", false, "none");
		int dms = 0;
		if(inputType.equals("bam")){
			dms=500;
		}else if(inputType.equals("xmap")){
			dms=10;
		}else{
			System.err.println("#<ERROR> inputType error\n");
		}
		minAln=tmp.getInt("minAln", false, dms);
		min=tmp.getInt("min", false, 500);//change from 2000 to 500
		max=tmp.getInt("max", false, 1000000);
		gq=5;//genotype quality
		FileWriter writer=new FileWriter(output);
		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
		LocalDateTime now = LocalDateTime.now();
		String infoPrint="";
		for(int i=0;i<args.length;i++){
			infoPrint+=args[i]+" ";
		}
		System.out.println("#<INFO-1> Arguments= "+infoPrint);
		System.out.println("#<INFO-2> Program starts at "+dtf.format(now));
		/////////////////do every thing in every chr
		final SamReaderFactory factory =
		         SamReaderFactory.makeDefault()
		             .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
		             SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		             .validationStringency(ValidationStringency.SILENT);

		final SamReader samReader = factory.open(new File(input));
		//samReader.getFileHeader().;
		SAMFileHeader samHead= samReader.getFileHeader();
		List<SAMSequenceRecord> samSequenceRecord = samHead.getSequenceDictionary().getSequences();
		Iterator<SAMSequenceRecord> iter2 = samSequenceRecord.iterator();
		ArrayList<String> allChr=new ArrayList<String>();
		while(iter2.hasNext()){
			allChr.add(iter2.next().getSequenceName());
		}
		samReader.close();
		//testModule();
		if(region.equals("all")){
			for(int i=0;i<allChr.size();i++){
				now = LocalDateTime.now();
				System.out.println("#<INFO-3> Running chromosome "+allChr.get(i)+" at "+dtf.format(now));
				runChr(allChr.get(i),shortBam,inputType,input,ir,args,writer);
			}
		}else{
			boolean testRegion=false;
			for(int i=0;i<allChr.size();i++){
				if(region.equals(allChr.get(i))){
					runChr(allChr.get(i),shortBam,inputType,input,ir,args,writer);
					testRegion=true;
				}
			}
			if(!testRegion){
				System.err.println("#<ERROR> input Region is not a chromosome in the samHeader!\n");
			}
		}
		writer.close();
		now = LocalDateTime.now();
		System.out.println("#<INFO-0> Program finishs at "+dtf.format(now));
	}
	private static void testModule() {
		// TODO Auto-generated method stub
		MergedSV tmp = null;
		for(int i=0;i<500;i++){
			for (int j=0;j<500;j++){
				int[] ref={i,0,0,0};
				int[] variant={j,0,0,0};
				double[] gen=tmp.genotype(ref, variant);
				double max;
				int maxID;
				double mid;
				int midID;
				if(gen[0]>gen[1]){
					if(gen[2]>gen[0]){
						max=gen[2];
						maxID=2;
						mid=gen[0];
						midID=0;
					}else{
						max=gen[0];
						maxID=0;
						if(gen[1]>gen[2]){
							mid=gen[1];
							midID=1;
						}else{
							mid=gen[2];
							midID=2;
						}
					}
				}else{
					if(gen[2]>gen[1]){
						max=gen[2];
						maxID=2;
						mid=gen[1];
						midID=1;
					}else{
						max=gen[1];
						maxID=1;
						if(gen[2]>gen[0]){
							mid=gen[2];
							midID=2;
						}else{
							mid=gen[0];
							midID=0;
						}
					}
				}
				int genotypeQuality=(int)(10*(max-mid));
				System.out.println(i+"\t"+j+"\t"+maxID+"\t"+genotypeQuality);
			}
		}
		
	}
	private static void runChr(String chr, String shortBam, String inputType, String input, String ir, String[] args, FileWriter writer) throws IOException {
		// TODO Auto-generated method stub	
		ArrayList<SV> sv = null;
		ReturnShortBam returnShortBam=null;
		ReturnLongBam returnLongBam=null;
		double[] errorAndStd=new double[6];
		if(!shortBam.equals("none")){
			DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
			LocalDateTime now = LocalDateTime.now();
			System.out.println("#<INFO-4> Progressing short "+shortBam + " for chromosome " + chr +" at "+dtf.format(now));
			returnShortBam=runShortBam(shortBam,chr);
		}
		if(inputType.equals("xmap")){
			sv=runXmap(input);
		}else if(inputType.equals("bam")){
			DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
			LocalDateTime now = LocalDateTime.now();
			System.out.println("#<INFO-4> Progressing chromosome " + chr +" at "+dtf.format(now));
			//outputErrorRate(input,chr);//just for testing
			errorAndStd=getAlignmentError(input,chr);
			returnLongBam=runBam(input,chr);
			
		}else if(inputType.equals("bamReadName")){
			//	sv=runBamRN(input);
		}else{
			System.err.println("#<ERROR> inputType Error="+inputType+"\n");
			System.exit(1);//error
		}
		//Not very good to put them here.
		ArrayList<MergedSV> mSV=Merge.MergeSignal(returnLongBam.sv,ir,chr);
//		mSV=tmp2.mergeComplex(mSV); //this function induces too many FDR
		mSV=RefDepth.removeLowCount(mSV);
		//update 1.1 and insert here.
		mSV=getDepth(mSV,input,chr,errorAndStd);
		//end update 1.1
//		mSV=getDepth(mSV,returnLongBam.depth);
		mSV=MergedSV.deletedRelate(mSV);
		if(!shortBam.equals("none")){
			ArrayList<MergedSV> clipSV=ClipRead.mSVwithClipRead(mSV,returnShortBam.signal,returnShortBam.depth);
			VCF.print(clipSV, args, writer,chr);	
		}
		VCF.print(mSV,args,writer,chr);// here sv to MergedSV		
	}
	//this function is to plot error rate
	private static void outputErrorRate(String input, String chr) {
		// TODO Auto-generated method stub
		final int start=107165550;
		final int end=107174550;
		final int interval=100;
		File tmpFile=new File(input);
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
		SAMRecordIterator iter=samReader.query(chr, start-100000, end+100000, true);
		SAMRecord sam;
		int[] refDepth=new int[4];
		double[] empty=new double[6];
		int id=1;
		while(iter.hasNext()){
			sam=iter.next();
			if(sam.getFlags()>=256){
				continue;
			}
			int startSite=((sam.getAlignmentStart()-50)/100+1)*100+50;
			int endSite=((sam.getAlignmentEnd()-50)/100)*100+50;
			for(int i=startSite;i<endSite;i+=100){
				int j=i+99;
				System.out.print(id+"\t"+i+"\t"+j+"\t");
				checkErrorRateFull(sam,i,i+99,empty);
			}
			id++;
		}

	}
	private static ArrayList<MergedSV> getDepth(ArrayList<MergedSV> mSV, String input, String chr,
			double[] errorAndStd) {
		// TODO Auto-generated method stub
		File tmpFile=new File(input);
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

		Iterator iSV=mSV.iterator();
		while(iSV.hasNext()){
			MergedSV msv=(MergedSV) iSV.next();
			SAMRecordIterator iter=samReader.query(msv.chr, msv.leftStart, msv.rightEnd, false);
			SAMRecord sam;
			int[] refDepth=new int[4];
			refDepth[0]=0;
			refDepth[1]=0;
			refDepth[2]=0;
			refDepth[3]=0;
			while(iter.hasNext()){
				sam=iter.next();
				if(sam.getFlags()>=256){
					continue;
				}
				if(sam.getAlignmentStart()<msv.leftStart && sam.getAlignmentEnd()>msv.rightEnd){
					//if(checkErrorRateFull(sam,msv.leftStart,msv.rightEnd,errorAndStd)){
					if(checkErrorRateFull(sam,msv.leftEnd,msv.rightStart,errorAndStd)){
						if(sam.getReadNegativeStrandFlag()){
							refDepth[1]++;
							refDepth[3]++;
						}else{
							refDepth[0]++;
							refDepth[2]++;
						}
					}
				}else if(sam.getAlignmentStart()<msv.leftStart && sam.getAlignmentEnd()>msv.leftEnd && sam.getAlignmentEnd()<=msv.rightEnd){ //if not full check left and right
					if(checkErrorRateLeft(sam,msv.leftEnd,msv.rightStart,errorAndStd)){
						if(sam.getReadNegativeStrandFlag()){
							refDepth[1]++;
						}else{
							refDepth[0]++;
						}
					}
				}else if(sam.getAlignmentStart()>=msv.leftStart && sam.getAlignmentStart()<msv.rightStart && sam.getAlignmentEnd()>msv.rightEnd){
					if(checkErrorRateRight(sam,msv.leftEnd,msv.rightStart,errorAndStd)){
						if(sam.getReadNegativeStrandFlag()){
							refDepth[3]++;
						}else{
							refDepth[2]++;
						}
					}
				}
			}
			iter.close();
			msv.setReference(refDepth);
		}
		return mSV;
	}
	private static boolean checkErrorRateRight(SAMRecord sam, int rightStart, int rightEnd, double[] errorAndStd) {
		// TODO Auto-generated method stub
		//System.out.print("right\t");
		String SAinfo=getSamMD(sam);		
		if(checkErrorRateRegion(rightStart,rightEnd,SAinfo,sam.getCigar(),sam.getAlignmentStart(),errorAndStd)){
			return true;
		}
		return false;
	}
	private static boolean checkErrorRateLeft(SAMRecord sam, int leftStart, int leftEnd, double[] errorAndStd) {
		// TODO Auto-generated method stub
		//System.out.print("left\t");
		String SAinfo=getSamMD(sam);		
		if(checkErrorRateRegion(leftStart,leftEnd,SAinfo,sam.getCigar(),sam.getAlignmentStart(),errorAndStd)){
			return true;
		}
		return false;
	}
	private static boolean checkErrorRateFull(SAMRecord sam, int start, int end, double[] errorAndStd) {
		// TODO Auto-generated method stub
		//System.out.print("full\t");
		String SAinfo=getSamMD(sam);		
		if(checkErrorRateRegion(start,end,SAinfo,sam.getCigar(),sam.getAlignmentStart(), errorAndStd)){
			return true;
		}
		return false;
	}
	private static String getSamMD(SAMRecord sam) {
		//System.out.print(sam.getAlignmentStart()+"\t"+sam.getAlignmentEnd()+"\t");

		//int readStart=sam.getReadPositionAtReferencePosition(start);
		//int readEnd  =sam.getReadPositionAtReferencePosition(end);		
		String line=sam.getSAMString();
		String[] lineSplit=line.split("\t");
		String SAinfo="";			
		for (int i=11;i<lineSplit.length;i++){ // start from 11 , skip other information
					// only check for SA:Z
					// should be updated for other method
			if(lineSplit[i].contains("MD:Z:")){
				SAinfo=lineSplit[i];
				return SAinfo;
			}
		}
		return SAinfo;
	}
	//version 1.1
	//key function to check Error rate
	//Insertion rate and deletion is got from cigar
	//Substitution rate is got from MD:Z it could only reconstruct reference but not the read because it does not record insertion
	private static boolean checkErrorRateRegion(int start, int end, String sAinfo, Cigar cigar,int aStart,double[] errorAndStd) {
		//sAinfo="MD:Z:4^A7^C4T6";
		String[] mdSplit=sAinfo.split("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");
//		System.out.println(Arrays.toString(mdSplit));
		int size=aStart;
		int i2=1;//new counter have to save
		double substitution=0;
		for(;i2<mdSplit.length;i2++){
			if(mdSplit[i2].contains("^")){
				size+=mdSplit[i2].length()-1;
			}else {
				try {
				    size+=Integer.parseInt(mdSplit[i2]);
				} catch (NumberFormatException e) {
				    size+=mdSplit[i2].length();
				}
			}
			if(size>start){
				break;
			}
		}
		for(;i2<mdSplit.length;i2++){
			if(mdSplit[i2].contains("^")){
				size+=mdSplit[i2].length()-1;
			}else {
				try {
				    size+=Integer.parseInt(mdSplit[i2]);
				} catch (NumberFormatException e) {
				    size+=mdSplit[i2].length();
				    substitution+=mdSplit[i2].length();
				}
			}
			if(size>=end){
				break;
			}
		}
		double snp=substitution/(end-start+1);
		//System.out.println("size="+size+"\tsnp="+snp);
		List<CigarElement> cigarList=cigar.getCigarElements();
		int ref=aStart;
		int read=0;
		double deletion=0;
		double insertion=0;
		for(int i=0;i<cigarList.size();i++){
			//D N ref++
			//I S read++
			//M = X ref++ read++
			//H P nothing++
			if(cigarList.get(i).getOperator()==CigarOperator.D || cigarList.get(i).getOperator()==CigarOperator.N){
				ref+=cigarList.get(i).getLength();
				if(ref>start && ref<=end){
					deletion+=cigarList.get(i).getLength();
				}
			}else if(cigarList.get(i).getOperator()==CigarOperator.I || cigarList.get(i).getOperator()==CigarOperator.S){
				read+=cigarList.get(i).getLength();
				if(ref>start && ref<=end && cigarList.get(i).getOperator()==CigarOperator.I){
					insertion+=cigarList.get(i).getLength();
				}
			}else if(cigarList.get(i).getOperator()==CigarOperator.M || cigarList.get(i).getOperator()==CigarOperator.X || cigarList.get(i).getOperator()==CigarOperator.EQ){
				read+=cigarList.get(i).getLength();
				ref+=cigarList.get(i).getLength();
			}else{
				//do nothing
			}
		}
		deletion=deletion/(end-start+1);
		insertion=insertion/(end-start+1);
		//System.out.println(snp+"\t"+deletion+"\t"+insertion);
		//System.out.println(errorAndStd[0]+"\t"+errorAndStd[1]+"\t"+errorAndStd[2]);
		//System.out.println("ref="+ref+"\tread="+read);
		if(snp>errorAndStd[0]+errorAndStd[1] || deletion>errorAndStd[2]+errorAndStd[3] || insertion>errorAndStd[4]+errorAndStd[5]){
		//if(snp>errorAndStd[0] || deletion>errorAndStd[2] || insertion>errorAndStd[4]){		
			return false;
		}else{
			return true;
		}
	}
	public static Integer tryParse(String text) {
		  try {
		    return Integer.parseInt(text);
		  } catch (NumberFormatException e) {
		    return text.length();
		  }
	}
	//input bam file and chr 
	//output error rate [0-2] snp. insertion ,deletion mean rate [3-5] snp, in-del std
	private static double[] getAlignmentError(String input, String chr) {
		// TODO Auto-generated method stub
		File tmpFile=new File(input);
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
		iter=samReader.query(chr, 1, 249250621, true);
		SAMRecord sam;
		int counter=0;
		int counterMax=100000;
		List<Double> snp = new ArrayList<Double>();
		List<Double> insertion = new ArrayList<Double>();
		List<Double> deletion = new ArrayList<Double>();
		while (iter.hasNext()){
			sam=iter.next();
			if(sam.getFlags()<256){
				counter++;
			}else{
				continue;
			}
			//you have to manually chop MD:Z to get this error rate.
			String line=sam.getSAMString();
			String[] lineSplit=line.split("\t");
			String SAinfo="";			
			for (int i=11;i<lineSplit.length;i++){ // start from 11 , skip other information
						// only check for SA:Z
						// should be updated for other method
				if(lineSplit[i].contains("NM:i:")){
					SAinfo=lineSplit[i];
					break;
				}
			}
			snp.add( Double.parseDouble(SAinfo.substring(5))/sam.getReadLength());
			Cigar cigar=sam.getCigar();
			List<CigarElement> cigarList=cigar.getCigarElements();
			double delete=0;
			double insert=0;
			for(int i=0;i<cigarList.size();i++){
				if(cigarList.get(i).getOperator()==CigarOperator.D){
					delete+=cigarList.get(i).getLength();
				}else if(cigarList.get(i).getOperator()==CigarOperator.I){
					insert+=cigarList.get(i).getLength();
				}
			}
			deletion.add(delete/sam.getReadLength());
			insertion.add(insert/sam.getReadLength());
			double s1=Double.parseDouble(SAinfo.substring(5))/sam.getReadLength();
			double d1=delete/sam.getReadLength();
			double i1=insert/sam.getReadLength();
			//System.out.println(sam.getAlignmentStart()+"\t"+sam.getAlignmentEnd()+"\t"+s1+"\t"+d1+"\t"+i1);
			if(counter>counterMax){
				break;
			}
		}
		double snpAll=0;
		double snpAll2=0;
		double insAll=0;
		double insAll2=0;
		double delAll=0;
		double delAll2=0;
		for(int i=0;i<snp.size();i++){
			snpAll+=snp.get(i);
			snpAll2+=snp.get(i)*snp.get(i);
			insAll+=insertion.get(i);
			insAll2+=insertion.get(i)*insertion.get(i);
			delAll+=deletion.get(i);
			delAll2+=deletion.get(i)*deletion.get(i);
			//System.out.println(snp.get(i)+"\t"+insertion.get(i)+"\t"+deletion.get(i));
		}
		double[] out=new double[6];
		out[0]=snpAll/snp.size();
		out[1]=Math.sqrt((snp.size()*snpAll2-snpAll*snpAll)/(snp.size()*snp.size()));

		out[2]=insAll/snp.size();
		out[3]=Math.sqrt((snp.size()*insAll2-insAll*insAll)/(snp.size()*snp.size()));

		out[4]=delAll/snp.size();
		out[5]=Math.sqrt((snp.size()*delAll2-delAll*delAll)/(snp.size()*snp.size()));
//		System.out.println(Arrays.toString(out));
//		double[] out2=new double[3];
//		for(int i=0;i<snp.size();i++){
//			out2[0]+=(snp.get(i)-out[0])*(snp.get(i)-out[0]);
//			out2[1]+=(insertion.get(i)-out[2])*(insertion.get(i)-out[2]);
//			out2[2]+=(deletion.get(i)-out[4])*(deletion.get(i)-out[4]);			
//		}
//		out2[0]=Math.sqrt(out2[0]/snp.size());
//		out2[1]=Math.sqrt(out2[1]/snp.size());
//		out2[2]=Math.sqrt(out2[2]/snp.size());
//		System.out.println(Arrays.toString(out2));
		return out;
	}
	private static void msvPrint(ArrayList<MergedSV> mSV) {
		// TODO Auto-generated method stub
		for(int i=0;i<mSV.size();i++){
			System.out.print(mSV.get(i).toString(i));
		}
	}
	private static void svPrint(ArrayList<SV> SV) {
		// TODO Auto-generated method stub
		for(int i=0;i<SV.size();i++){
			System.out.print(SV.get(i).toString());
		}
	}
	/*
	private static ArrayList<SV> runBamRN(String input) {
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.LENIENT);
		SAMFileReader samReader = new SAMFileReader(new File(input));	
		String line;
		SAMRecordIterator iter=samReader.iterator();
		SAMRecord sam;
		///here is the start
		sam=iter.next();
		String readName=sam.getReadName();
		if(!depth.containsKey(sam.getReferenceName())){
			depth.put(sam.getReferenceName(), new ArrayList<startEnd>());
		}
		depth.get(sam.getReferenceName()).add(new startEnd(sam.getAlignmentStart(),sam.getAlignmentEnd(),sam.getReadNegativeStrandFlag()));
		String strand;
		if(sam.getReadNegativeStrandFlag()){
			strand="-"; // checked
		}else{
			strand="+"; // checked
		}
		ArrayList<Alignment> all=new ArrayList<Alignment>();
		all.add(new Alignment("bam",sam.getReferenceName(),sam.getAlignmentStart(),sam.getCigarString(),strand));	
		//loop start
		Read read;
		while (iter.hasNext()){
			sam=iter.next();
			line=sam.getSAMString();
			if(!depth.containsKey(sam.getReferenceName())){
				depth.put(sam.getReferenceName(), new ArrayList<startEnd>());
			}
			depth.get(sam.getReferenceName()).add(new startEnd(sam.getAlignmentStart(),sam.getAlignmentEnd(),sam.getReadNegativeStrandFlag()));
			if(sam.getReadNegativeStrandFlag()){
				strand="-"; // checked
			}else{
				strand="+"; // checked
			}
			if(sam.getReadName().equals(readName)){
				all.add(new Alignment("bam",sam.getReferenceName(),sam.getAlignmentStart(),sam.getCigarString(),strand));	
			}else{
				read=new Read(all,minAln);// filter Alignment here
				if(read.getAlignmentArray().size()>1){
					for(int i=1;i<read.getAlignmentArray().size();i++){
						SV tmp=new SV(read.getAlignmentArray().get(i-1),read.getAlignmentArray().get(i));
						if(tmp.getSVcorrect()){
							sv.add(tmp);
						}
					}
				}
				readName=sam.getReadName();
				all.clear();
				all.add(new Alignment("bam",sam.getReferenceName(),sam.getAlignmentStart(),sam.getCigarString(),strand));	
			//end of the reading bam
			}
		}
		read=new Read(all,minAln);
		if(read.getAlignmentArray().size()>1){
			for(int i=1;i<read.getAlignmentArray().size();i++){
				SV tmp=new SV(read.getAlignmentArray().get(i-1),read.getAlignmentArray().get(i));
				if(tmp.getSVcorrect()){
					sv.add(tmp);
				}
			}
		}
		samReader.close();
		Collections.sort(sv);
		return sv;
	}
	*/
	private static ArrayList<MergedSV> getDepth(ArrayList<MergedSV> mSV, startEnd[] depth) {
		//for all SV
		for(int i=0;i<mSV.size();i++){
			//for this SV to all depth in this chromosome
			getDepthSingle(mSV.get(i),depth);// actually mSV is already updated
		}
		return mSV;
	}
	public static void getDepthSingle(MergedSV MergedSV, startEnd[] depth) {
		// TODO Auto-generated method stub
		int[] ref= new int[4];
		int[] ref1=getReadCount(MergedSV.getLeftStart(),MergedSV.getLeftEnd(),depth);
		//depth.get(chr) is sorted!!
		int[] ref2=getReadCount(MergedSV.getRightStart(),MergedSV.getRightEnd(),depth);
		ref[0]=ref1[0];
		ref[1]=ref1[1];
		ref[2]=ref2[0];
		ref[3]=ref2[0];
		MergedSV.setReference(ref);
	}
	private static int[] getReadCount(int start, int end, startEnd[] depth) {
		// TODO Auto-generated method stub
		int[] result=new int[2];
		startEnd tmp=null;
		for(int j=0;j<depth.length;j++){
			tmp=depth[j];
			if(tmp.start<start && end<tmp.end){
				if(tmp.negativeStrandFlag){
					result[1]++;
				}else{
					result[0]++;
				}
			}else if(end<tmp.start){
				break;
			}
		}
		return result;
	}
	private static ReturnShortBam runShortBam(String input,String chr) {
		ArrayList<startEnd> depthShortList=new ArrayList<startEnd>();
		ArrayList<ClipRead> clipList=new ArrayList<ClipRead>();
		final SamReaderFactory factory =
		         SamReaderFactory.makeDefault()
	             .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
	             SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	             .validationStringency(ValidationStringency.SILENT);

		final SamReader samReader = factory.open(new File(input));
		SAMRecordIterator iter;
		iter=samReader.query(chr, 1, 249250621, true);
		int runCounter=0;
		String cigar;
		String[] digital;
		String[] word;
		SAMRecord sam;
		while (iter.hasNext()){
			sam=iter.next();
			if(sam.getFlags()<256){
				depthShortList.add(new startEnd(sam.getAlignmentStart(),sam.getAlignmentEnd(),sam.getReadNegativeStrandFlag()));
			}			
			cigar=sam.getCigarString();
			digital=cigar.split("M|I|D|N|S|H|P|=|X");
			word=cigar.split("\\d+");
			//System.out.println(cigar+"sizeOfWord"+word.length+"\tword="+Arrays.toString(word));
			//word start from 1 word[0]="";
			if(word.length>=2 && word[1].equals("S") && Integer.parseInt(digital[0])>clipSize){
				clipList.add(new ClipRead(sam.getAlignmentStart(),true));
			}
			if(word[word.length-1].equals("S") && Integer.parseInt(digital[digital.length-1])>clipSize && sam.getAlignmentEnd()!=0){
				clipList.add(new ClipRead(sam.getAlignmentEnd(),false));
			}
			
			//end of the reading bam
			runCounter++;
			if((runCounter % 50000000)==0){
//				DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
//				LocalDateTime now = LocalDateTime.now();
				System.gc();//manually release memory !!!!!
//				System.out.println(runCounter/1000000+"M reads have been read at "+dtf.format(now)+" with memory "+Runtime.getRuntime().totalMemory()+" for ref "+sam.getReferenceName());
			}
			sam.clearAttributes();
		}
		try {
			samReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		startEnd[] depthShort = list2array(depthShortList);
		ClipRead[] clipRaw=new ClipRead[clipList.size()];
		clipRaw=clipList.toArray(clipRaw);
		ClipRead[] clip = ClipRead.mergeAll(clipRaw);
		ReturnShortBam result=new ReturnShortBam(depthShort,clip,chr);
		return result;
	}
		
	private static ReturnLongBam runBam(String input, String chr) throws IOException {		
		final SamReaderFactory factory =
	         SamReaderFactory.makeDefault()
             .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
             SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
             .validationStringency(ValidationStringency.SILENT);

		final SamReader samReader = factory.open(new File(input));
		String line;
		SAMRecordIterator iter;
		iter=samReader.query(chr, 1, 2147483647, true);//max int for this chromosome
		//Not using it anymore
		ArrayList<startEnd> depthList=new ArrayList<startEnd>();
		//
		ArrayList<SV> sv=new ArrayList<SV>();
		SAMRecord sam;
		// too much memory
		//HashMap<String, Integer> mul=new HashMap<String, Integer>();//just make sure they are analysis
		while (iter.hasNext()){
			sam=iter.next();
			Cigar c=sam.getCigar();
			line=sam.getSAMString();
			if(sam.getFlags()<256){
				startEnd tmp=new startEnd(sam.getAlignmentStart(),sam.getAlignmentEnd(),sam.getReadNegativeStrandFlag());
				depthList.add(tmp);
			}
			if(line.contains("SA:Z")){
				if(sam.getFlags()>=256){
					//do nothing
				}else{
					//mul.put(sam.getReadName(), 1);
					String[] lineSplit=line.split("\t");
					String SAinfo="";			
					for (int i=11;i<lineSplit.length;i++){ // start from 11 , skip other information
								// only check for SA:Z
								// should be updated for other method
						if(lineSplit[i].contains("SA:Z")){
							SAinfo=lineSplit[i];
							break;
						}
					}
					//System.out.println(SAinfo);
					String strand;
					if(sam.getReadNegativeStrandFlag()){
						strand="-"; // checked
					}else{
						strand="+"; // checked
					}
					//System.out.println(sam.getSAMString()+strand+" "+sam.getCigarString());
					String[] mulAlign=SAinfo.split("SA:Z:|;");
					ArrayList<Alignment> all=new ArrayList<Alignment>();
					all.add(new Alignment("bam",sam.getReferenceName(),sam.getAlignmentStart(),sam.getCigarString(),strand));					
					String[] mulInfo;
					for(int i=1;i<mulAlign.length;i++){
						mulInfo=mulAlign[i].split(",");
						all.add(new Alignment("bam",mulInfo[0],Integer.parseInt(mulInfo[1]),mulInfo[3],mulInfo[2]));
					}
					Read read=new Read(all,minAln);// filter Alignment here
					if(read.getAlignmentArray().size()>1){
						for(int i=1;i<read.getAlignmentArray().size();i++){
							SV tmp=new SV(read.getAlignmentArray().get(i-1),read.getAlignmentArray().get(i));
							if(tmp.getSVcorrect()){
								sv.add(tmp);
						//		System.out.print(tmp.toString(10086));
							}else{
						//		System.out.println("false for"+read.getAlignmentArray().get(i-1).getMappingStart()+"\t"+read.getAlignmentArray().get(i).getMappingStart());
							}
						}
					}
					if(test){
						System.out.println("end of this read");
					}
				}
				//end of if
			}
			sam.clearAttributes();
			//end of the reading bam
		}
		startEnd[] depth= new startEnd[depthList.size()];
		depth=list2array(depthList);
		samReader.close();
		Collections.sort(sv);
		if(test){
			Iterator svIter=sv.iterator();
			while(svIter.hasNext()){
				System.out.print(svIter.next().toString());
			}
		}
		// Do something for depth and sv to get the mapping error
		
		//
		ReturnLongBam returnLongBam=new ReturnLongBam(depth,sv,chr);
		return returnLongBam;
	}
	//convert list to array to save time
	private static startEnd[] list2array(ArrayList<startEnd> depthList) {
		startEnd[] result =new startEnd[depthList.size()];
		result=depthList.toArray(result);
		return result;
	}
	private static ArrayList<SV> runXmap(String input) throws FileNotFoundException {
		// TODO Auto-generated method stub
		BufferedReader inputReader = new BufferedReader(new FileReader(input));
		String line;
		HashMap<String, String> all= new HashMap<String, String>(); // hash <contigID,line1\tline2\t..lineN\t>
		HashMap<String, Integer> count= new HashMap<String, Integer>(); // hash <contigID,numberOfLine>
		String[] info;
	    try {
			while ((line = inputReader.readLine()) != null) {
				if(! line.subSequence(0, 1).equals("#")){
					info=line.split("\t");
					if(all.containsKey(info[1])){
						all.put(info[1], all.get(info[1])+line+"\n");
						count.put(info[1], count.get(info[1])+1);
					}else{
						all.put(info[1] , line+"\n");
						count.put(info[1],1);
					}
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    Enumeration<String> keys=(Enumeration<String>) all.keySet();
	    String[] lineInfo;
	    String[] oneLineInfo;
	    String[] numberInfo;
	    while (keys.hasMoreElements()) {
	        String key = (String) keys.nextElement();
	       if(true){
	        //if(count.get(key)>1){ for n=1
	        		lineInfo=all.get(key).split("\n");
	        		int[][] query=new int[lineInfo.length][];
	        		for (int i=0;i<lineInfo.length;i++){
//////////////////////////////////////////////////////////////////////////////////////////
//////////IF YOU WANT TO USE READ MODIFY here///////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
	        			oneLineInfo=lineInfo[i].split("\t");
	        			numberInfo=oneLineInfo[13].split("\\(|\\,|\\)");// in here simply split
	        			// (ref1,query1)(ref2,query2)...(refN,queryN) 
	        			//after split==> ref1=x[1] query1=x[2] ref2=4 query=5 ... 
	        			//0 3 are gaps...
	        			query[i]=new int[numberInfo.length/3];
	        			for (int j=0;j*3<numberInfo.length;j++){
	        				if(j*3+3<numberInfo.length){
	        					if(Math.abs(Integer.parseInt(numberInfo[j*3+1])-Integer.parseInt(numberInfo[j*3+4]))>minAln 
	        							&& Math.abs(Integer.parseInt(numberInfo[j*3+2])-Integer.parseInt(numberInfo[j*3+5]))>minAln){
	        						System.out.println(numberInfo[j*3+1]+"~"+numberInfo[j*3+4]+"--"+Integer.parseInt(numberInfo[j*3+2])+"~"+Integer.parseInt(numberInfo[j*3+5]));
	        						System.out.println(lineInfo[i]);
	        					}
	        				}
	        				query[i][j]=Integer.parseInt(numberInfo[j*3+1]);
	        			};
	        			//System.out.println(Arrays.toString(numberInfo));
	        		}
	        		int[][] startEnd=new int[lineInfo.length][2];
	        		for (int i=0;i<lineInfo.length;i++){
	        			startEnd[i][0]=query[i][0];
	        			startEnd[i][1]=query[i][query[i].length-1];
	        		}
	        		ArrayList<Integer> remainAlignment=minimumOverlap(startEnd);
	        		if(remainAlignment.size()>1){
	        			int value;
	        			for(int q=0;q<remainAlignment.size();q++){
	        				value=remainAlignment.get(q);
	        				//System.out.printf("q=%d\tstart=%d\tend=%d\n",startEnd[remainAlignment.get(q)][0],startEnd[remainAlignment.get(q)][1]);
	        				//System.out.println(lineInfo[value]);
	        			}
	        			//System.out.println();
	        		//System.out.println(Arrays.toString(remainAlignment.toArray()));
	        		}
	        	}
	    }
	    ArrayList<SV> sv = null;
	    return sv;
	}

	private static ArrayList<Integer> minimumOverlap(int[][] startEnd) {
		// time O(n**2)
		ArrayList<Integer> remainAlignment = new ArrayList<Integer>();
		boolean keep;
		for (int i=0;i<startEnd.length;i++){
			keep=true;
			for (int j=0;j<startEnd.length;j++){
				if(j!=i){
					//System.out.println(startEnd[j][0]+" "+startEnd[i][0]+" "+startEnd[i][1]+" "+startEnd[j][1]);
					if((startEnd[j][0]<=startEnd[i][0] && startEnd[i][1]<startEnd[j][1] )||
							(startEnd[j][0]<startEnd[i][0] && startEnd[i][1]<=startEnd[j][1])){ // double largest ???
						keep=false;
						break;
					}else{
					}
				}
			}
			if (keep){
				remainAlignment.add(i);
			}
		}
		// check whether there is double largest value like (2,10) (2,10) (6,90), remove duplicated one(keep the last one).
		for(int i=0;i<remainAlignment.size();i++){
			for(int j=i+1;j<remainAlignment.size();j++){
				if(startEnd[remainAlignment.get(i)][0]==startEnd[remainAlignment.get(j)][0] 
						&& startEnd[remainAlignment.get(i)][1]==startEnd[remainAlignment.get(j)][1]){
					remainAlignment.remove(i);
				}
			}
		}
		// TODO Auto-generated method stub
		return remainAlignment;
	}

}
