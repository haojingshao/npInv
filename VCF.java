package inversion;


import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
//import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Iterator;

public class VCF {
	public static double version=Main.version;
	public static boolean headerPrintOrNot=false;
	public static String header="";
	public static final int mergedSize=Main.mergedSize;
	public static int id=1;
	public static void print(ArrayList<MergedSV> svList,String[] args, FileWriter writer, String chr)throws IOException{	
		if(!headerPrintOrNot){
			writer.write(header);
			headerPrintOrNot=true;
		}
		for(int i=0;i<svList.size();i++){
			if(svList.get(i).chr.equals(chr)){ // delete other chromosome inversion from this chromosome analysis
				writer.write(svList.get(i).toString(id));// start from 1 for ID
				id++;
			}
		}
	}
	static void createHeader(String[] args, ArrayList<String> allChr) {
//		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
//		LocalDateTime now = LocalDateTime.now();
		header+="##fileformat=VCFv4.2\n";
//		header+="##fileDate="+dtf.format(now)+"\n";
		header+="##arguments="+String.join(" ",args)+"\n";
		header+="##method=npInv"+Main.version+"\n";
		for(int i=0;i<allChr.size();i++){
			header+="##contig=<ID="+allChr.get(i)+">\n";
		}
		header+="##ALT=<ID=INV,Description=\"Inversion\">\n";
		//header+="##ALT=<ID=MMEJ,Description=\"Micro-homology end joining\">\n";
		header+="##INFO=<ID=Ref,Description=\"Genotype is reference\">\n";
		header+="##INFO=<ID=Single,Description=\"Inversion supporting read only from single side or from single strand\">\n";
		header+="##INFO=<ID=LowGQ,Description=\"Genotype quality < 20\">\n";
		header+="##INFO=<ID=Long,Description=\"Inversion size is longer than expected(Default:1M)\">\n";
		header+="##INFO=<ID=Short,Description=\"Inversion size is shorter than expected(Default:2k)\">\n";
		header+="##INFO=<ID=PASS,Description=\"High quality result\">\n";	
		header+="##INFO=<ID=END,Number=1,Type=Integer,Description=\"Average end position of the inversion\">\n";
		header+="##INFO=<ID=IR,Number=1,Type=Integer,Description=\"Size of inverted repeat\">\n";
		header+="##INFO=<ID=MH,Number=1,Type=Integer,Description=\"Size of micro-homology\">\n";
		header+="##INFO=<ID=MECHANISM,Number=1,Type=String,Description=\"NHEJ=Non homology end joining, NAHR=Non allelic homology recombination (IR size>500bp)\">\n";
		//header+="##INFO=<ID=NAHR,Description=\"Non allelic homology recombination (IR size>500bp)\">\n";
		header+="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		header+="##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
		header+="##FORMAT=<ID=DR,Number=4,Type=Integer,Description=\"Depth reference for "
				+ "left breakpoint forward, left breakpoint reverse, right breakpoint forward and right breakpoint reverse\">\n";
		header+="##FORMAT=<ID=DV,Number=4,Type=Integer,Description=\"Depth variant for "
				+ "left breakpoint forward, left breakpoint reverse, right breakpoint forward and right breakpoint reverse\">\n";
		header+="##FORMAT=<ID=BP,Number=4,Type=Integer,Description=\"BreakPoint for minimum left start,maximum left end,minimum right start and maximum right end\">\n";
		header+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
	}
}
