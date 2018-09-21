package inversion;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;
import java.lang.*;

public class SV implements Comparable<SV> {
	public final static boolean test=Main.test;
	public String type;
	public String chr;
	public String chr2;// the same for deletion insertion 
	public int start1;
	public int start2;// start2 is the end for deletion
//	public int svSize;
	public boolean correct;// true=right false=correct
	//public String orientation;
	public int orientation;
	public int homo;
	// May be you should not use the constructor to check it is a SV or not
	public SV(Alignment first, Alignment second){
		this.correct=false;
		boolean firstStrand=first.getStrand().equals("+");
		boolean secondStrand=second.getStrand().equals("+");
		if(first.getChr().equals(second.getChr())){
			this.chr=first.getChr();
			this.chr2=second.getChr();
			if((first.getStrand().equals("+") && second.getStrand().equals("+") ) ||
					(first.getStrand().equals("-") && second.getStrand().equals("-") )){
			}else if((first.getStrand().equals("+") && second.getStrand().equals("-") ) ||
					(first.getStrand().equals("-") && second.getStrand().equals("+") )){
				// different orientation
				this.type="INV";
				if(Math.abs(first.getRefStart()-second.getRefStart())<100 || Math.abs(first.getRefEnd()-second.getRefEnd())<100){
	//				System.out.println("false abs");
					this.correct=false;// reverse read
				}
				int tmp=first.getMappingEnd()-second.getMappingStart();
				if(tmp<0){
					this.homo=0;
				}else {
					this.homo=tmp;
				}
				//this is choosing outer of inversion ,not to deleting the IR
				//4type
				//Draw picture
				if(firstStrand){
					if(first.getRefEnd()<second.getRefEnd()){
						//type 0: 
						this.start1=first.getRefEnd()-this.homo;
						this.start2=second.getRefEnd();
						this.orientation=0;
						//this.orientation="LF";
					}else{
						this.start1=first.getRefEnd();
						this.start2=second.getRefEnd()-this.homo;
						this.orientation=1;
						//this.orientation="LR";
					}
				}else if(secondStrand){
					if(first.getRefStart()<second.getRefStart()){
						this.start1=first.getRefStart();
						this.start2=second.getRefStart()+this.homo;
						this.orientation=2;
						//this.orientation="RF";
					}else{
						this.start1=first.getRefStart()+this.homo;
						this.start2=second.getRefStart();
						this.orientation=3;
						//this.orientation="RR";
					}
				}else{
					System.err.println("Program error for firstStrand/secondStrand\n");
					System.exit(1);
				}
				this.correct=true;
			}else{
				System.err.println("correct in orientation."+first.getStrand()+"\t"+ second.getStrand()+"Contact author\n");
				System.exit(1);
			}
		}else{
		}
		// sort the chromosome
		if(first.getChr().equals(second.getChr())){
			if(this.getStart1()>this.getStart2()){
				int tmp=this.getStart1();
				this.start1=this.getStart2();
				this.start2=tmp;
			}
		}
	}
	public SV(String string, String chr, int start, int end, int size) {
		// TODO Auto-generated constructor stub
		// not use
		this.type=string;
		this.chr=chr;
		this.chr2=chr;
		this.start1=start;
		this.start2=end;
//		this.svSize=size;
		this.correct=true;//right
//		System.out.println(this.toString());
	}
	public SV(String string) {
		//this.sv+="DEL\t"+this.chr+"\t"+ Integer.toString(start)+"\t"+ 
		//Integer.toString(end)+"\t"+ digital[i]+"\n";
		String[] info=string.split("\t");
		this.type=info[0];
		this.chr=info[1];
		this.chr2=info[2];
		this.start1=Integer.parseInt(info[3]);
		this.start2=Integer.parseInt(info[4]);
		//this.svSize=Integer.parseInt(info[5]);
	}
	public boolean getSVcorrect(){
		return this.correct;
	}

	@Override
	public int compareTo(SV o) {
		// TODO Auto-generated method stub
		int out=compareChr(this.chr,o.chr);
		if(out==0){
			return this.start1-o.start1;
		}else{
			return out;
		}
	}
	public static int compareChr(String a,String b){	
		return a.compareTo(b);
	}
//	public static int compareChr(String a,String b){		
//		String[] a1=a.split("(?=\\d)(?<!\\d)");
//		String[] b1=b.split("(?=\\d)(?<!\\d)");
//		if(!a1[0].equals(b1[0])){
//			return a.compareTo(b);
//		}else{
//			if(a1.length==2 && b1.length==2){
//				try{
//					return Integer.parseInt(a1[1])-Integer.parseInt(b1[1]);
//				}catch(NumberFormatException e){
//					return a.compareTo(b);
//				}
//			}else{
//				return a.compareTo(b);
//			}
//		}
//	}
	public String toString(){
		String info="CHR2="+this.chr2+";END="+this.start2+";HOMO="+this.homo+";orientation="+this.orientation+";";
		String out=this.chr+"\t"+this.start1+" \t"+0+"\tREF\t"+this.type+"\tQUAL\tFILTER\t"+info+"\tGT\tNA\n";
		return out;	
	}
	public String toString(int i) {
		// TODO Auto-generated method stub
		String info="CHR2="+this.chr2+";END="+this.start2+";HOMO="+this.homo+";orientation="+this.orientation+";";
		String out=this.chr+"\t"+this.start1+" \t"+i+"\tREF\t"+this.type+"\tQUAL\tFILTER\t"+info+"\tGT\tNA\n";
		return out;
	}
	public String getType(){
		return this.type;
	}
	public String getchr(){
		return this.chr;
	}
	public String getChr2(){
		return this.chr2;
	}
	public int getStart1(){
		return this.start1;
	}
	public int getStart2(){
		return this.start2;
	}
	public int getOrientation(){
		return this.orientation;
	}
	public int getHomo(){
		return this.homo;
	}
}
