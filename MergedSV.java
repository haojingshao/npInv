package inversion;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

public class MergedSV implements Comparable<MergedSV>{
	public boolean KeepOrNot=true;
	public String type;
	public boolean complexOrNot;
	public String chr;
	public int start1;
	public int start2;// start2 is the end for deletion
//	private int svSize;
	public int[] orientation=new int[4];
	public int[] reference=new int[4];
//	public static Hashtable<String, Integer> o2int= new Hashtable<String, Integer>(){
//		private static final long serialVersionUID = 1L;
//	{put("LF",0);put("LR",1);put("RF",2);put("RR",3);}};
	public static double variantRate=Main.variantRate;
	public int homo;
	public int homoMax=0;
	public int leftStart;
	public int leftEnd;
	public int rightStart;
	public int rightEnd;
	public int numberOfRead;
	public int genotype=0; //0="0/0" 1="0/1" 2="1/1"
	public int genotypeQuality=0;
	public MergedSV(ArrayList<SV> sv){
		this.complexOrNot=false;
		this.numberOfRead=sv.size();
		this.type=sv.get(0).getType();
		this.chr=sv.get(0).getchr();
		this.leftStart=sv.get(0).getStart1();
		this.leftEnd=sv.get(0).getStart1()+sv.get(0).getHomo();
		this.rightStart=sv.get(0).getStart2()-sv.get(0).getHomo();
		this.rightEnd=sv.get(0).getStart2();		
		int countHomo=0;
		double dstart1=0;
		double dstart2=0;
		for(int i=0;i<sv.size();i++){
			this.orientation[sv.get(i).orientation]+=1;
			dstart1+=sv.get(i).getStart1();
			dstart2+=sv.get(i).getStart2();
			if(sv.get(i).getHomo()>0){
				countHomo++;
				this.homo+=sv.get(i).getHomo();
				if(sv.get(i).getHomo()>this.homoMax){
					this.homoMax=sv.get(i).getHomo();
				}
			}
			if(sv.get(i).getStart1()<this.leftStart){
				this.leftStart=sv.get(i).getStart1();
			}
			if(sv.get(i).getStart1()+sv.get(i).getHomo()>this.leftEnd){
				this.leftEnd=sv.get(i).getStart1()+sv.get(i).getHomo();
			}
			if(sv.get(i).getStart2()-sv.get(i).getHomo()<this.rightStart){
				this.rightStart=sv.get(i).getStart2()-sv.get(i).getHomo();
			}
			if(sv.get(i).getStart2()>this.rightEnd){
				this.rightEnd=sv.get(i).getStart2();	
			}
		}
		this.start1=(int) (dstart1/sv.size());
		this.start2=(int) (dstart2/sv.size());
		if(countHomo>0){
			this.homo=this.homo/countHomo;
		}else{
			this.homo=0;
		}
		if(this.homo>=500){
			this.type="NAHR";
		}else{
			this.type="NHEJ";
		}
	}
	// for complexSV
	public MergedSV(MergedSV MergedSV, MergedSV MergedSV2) {
		this.complexOrNot=true;
		this.type="Complex";
		this.chr=MergedSV.chr;
		this.start1=(MergedSV.start1+MergedSV2.start1)/2;
		this.start2=(MergedSV.start2+MergedSV2.start2)/2;
		this.orientation[0]=MergedSV.orientation[0]+MergedSV2.orientation[0];
		this.orientation[1]=MergedSV.orientation[1]+MergedSV2.orientation[1];
		this.orientation[2]=MergedSV.orientation[2]+MergedSV2.orientation[2];
		this.orientation[3]=MergedSV.orientation[3]+MergedSV2.orientation[3];
		this.homo=(MergedSV.homo+MergedSV2.homo)/2;
		this.homoMax=(MergedSV.homoMax>MergedSV2.homoMax)?MergedSV.homoMax:MergedSV2.homoMax;
		this.leftStart=(MergedSV.leftStart<MergedSV2.leftStart)?MergedSV.leftStart:MergedSV2.leftStart;
		this.leftEnd=(MergedSV.leftEnd>MergedSV2.leftEnd)?MergedSV.leftEnd:MergedSV2.leftEnd;
		this.rightStart=(MergedSV.rightStart<MergedSV2.rightStart)?MergedSV.rightStart:MergedSV2.rightStart;
		this.rightEnd=(MergedSV.rightEnd>MergedSV2.rightEnd)?MergedSV.rightEnd:MergedSV2.rightEnd;
		this.numberOfRead=MergedSV.numberOfRead+MergedSV2.numberOfRead;
	}
	public MergedSV(ArrayList<SV> sv, IRPair ir) {
		// TODO Auto-generated constructor stub
		this.complexOrNot=false;
		this.type="NAHR-IR";
		this.chr=sv.get(0).chr;
		this.start1=(ir.leftEnd+ir.leftStart)/2;
		this.start2=(ir.rightEnd+ir.rightStart)/2;
		this.homo=(ir.leftEnd-ir.leftStart+1+ir.rightEnd-ir.rightStart+1)/2;
		this.homoMax=(ir.leftEnd-ir.leftStart>ir.rightEnd-ir.rightStart)?ir.leftEnd-ir.leftStart+1:ir.rightEnd-ir.rightStart+1;
		this.leftStart=ir.leftStart;// increase genotype error because the software will align more
		this.leftEnd=ir.leftEnd;
		this.rightStart=ir.rightStart;
		this.rightEnd=ir.rightEnd;
		MergedSV tmp=new MergedSV(sv);
		this.leftStart=(tmp.leftStart<ir.leftStart)?tmp.leftStart:ir.leftStart;
		this.leftEnd=(tmp.leftEnd>ir.leftEnd)?tmp.leftEnd:ir.leftEnd;
		this.rightStart=(tmp.rightStart<ir.rightStart)?tmp.rightStart:ir.rightStart;
		this.rightEnd=(tmp.rightEnd>ir.rightEnd)?tmp.rightEnd:ir.rightEnd;	
		this.numberOfRead=sv.size();
		for(int i=0;i<sv.size();i++){
			this.orientation[sv.get(i).orientation]++;
		}
	}
	//for clip reads
	public MergedSV(MergedSV tmpMSV, ArrayList<ClipRead> tmpClip1, ArrayList<ClipRead> tmpClip2) {
		// TODO Auto-generated constructor stub
		//they are aleady inside
		this.complexOrNot=false;
		this.type="Short";
		this.chr=tmpMSV.chr;
		this.homo=tmpMSV.homo;
		this.homoMax=tmpMSV.homoMax;
		this.start1=tmpMSV.start1;
		this.start2=tmpMSV.start2;
		boolean[] signal=new boolean[4];
		Arrays.fill(signal, false);
		int tf=0;
		int tr=0;
		int counterF=0;
		int counterR=0;
		for(int i=0;i<tmpClip1.size();i++){
			ClipRead tmp=tmpClip1.get(i);
			if(tmp.forwardOrNot){
				signal[0]=true;
				tf+=tmp.pos*tmp.count;
				counterF+=tmp.count;
			}else{
				signal[1]=true;
				tr+=tmp.pos*tmp.count;
				counterR+=tmp.count;
			}
		}
		this.leftStart=(counterF>0)?tf/counterF:0;
		this.leftEnd=(counterR>0)?tr/counterR:0;
		this.numberOfRead=counterF+counterR;
		this.orientation[0]=counterF;
		this.orientation[1]=counterR;
		tf=0;
		tr=0;
		counterF=0;
		counterR=0;
		for(int i=0;i<tmpClip2.size();i++){
			ClipRead tmp=tmpClip2.get(i);
			if(tmp.forwardOrNot){
				signal[2]=true;
				tf+=tmp.pos*tmp.count;
				counterF+=tmp.count;
			}else{
				signal[3]=true;
				tr+=tmp.pos*tmp.count;
				counterR+=tmp.count;
			}
		}
		this.rightStart=(counterF>0)?tf/counterF:0;
		this.rightEnd=(counterR>0)?tr/counterR:0;
		this.numberOfRead+=counterF+counterR;
		this.orientation[2]=counterF;
		this.orientation[3]=counterR;
		if(signal[0] && signal[1] && signal[2] && signal[3]){
			this.KeepOrNot=true;
		}else{
			this.KeepOrNot=false;
		}
	}
	public String toString(int i){
		String qual=".";
		String filter="LowQual";
		if(this.complexOrNot){
			this.type="Complex";
		}
		//if(this.genotypeQuality>=20 && (this.orientation[0]+this.orientation[1])>0 && (this.orientation[2]+this.orientation[3])>0){
		int count=0;
		if (this.orientation[0]>0){count++;}
		if (this.orientation[1]>0){count++;}
		if (this.orientation[2]>0){count++;}
		if (this.orientation[3]>0){count++;}
		//GQ>20 
		//inversion Count>3
		//genotype !=0/0
		//size<1M
		//size>2k
		if(this.genotype==0){
			filter="Ref";
		}else if(count<3){
			filter="Single";
		}else if(this.rightEnd - this.leftStart>=Main.max){
			filter="Long";
		}else if(this.start2 - this.start1<=Main.min){
			filter="Short";
		}else if(this.genotypeQuality<Main.gq){
			filter="LowGQ";
		}else if(this.genotypeQuality>=Main.gq && count>=3 && this.genotype!=0 && this.rightEnd - this.leftStart<Main.max && this.start2 - this.start1>Main.min){
			filter="PASS";
		}
		String GT="";
		if (this.genotype==0){
			GT="0/0";
		}else if(this.genotype==1){
			GT="0/1";
		}else if(this.genotype==2){
			GT="1/1";
		}else{
			GT="NA";//this is error!
		}
		String info="END="+this.start2;
		if(this.type.equals("NAHR") || this.type.equals("NAHR-IR")){
			info+=";IR="+this.homo;
			/////////////// for the headline
			info+=";MECHANISM="+this.type;
			this.type="INV";
		}else if(this.type.equals("NHEJ")){
			////
			info+=";MECHANISM="+this.type;
			this.type="INV";
		}else if(this.type.equals("MMEJ")){
			info+=";MH="+this.homo;
			////
			info+=";MECHANISM="+this.type;
			this.type="INV";
		}else if(this.complexOrNot){
			if(this.leftEnd-this.leftStart>this.rightEnd-this.rightStart){
				info+=";Deletion="+this.leftStart+"-"+this.leftEnd;
			}else{
				info+=";Deletion="+this.rightStart+"-"+this.rightEnd;
			}
		}
		//+";GT="+this.genotype+";QT="+this.genotypeQuality
		//		";SESE="+this.leftStart+"|"+this.leftEnd+"|"+this.rightStart+"|"+this.rightEnd+";";
		String format="GT:GQ:DR:DV:BP";
		String formatResult=GT+":"+this.genotypeQuality+":"
				+this.reference[0]+","+this.reference[1]+","+this.reference[2]+","+this.reference[3]
				+":"+this.orientation[0]+","+this.orientation[1]+","+this.orientation[2]+","+this.orientation[3]
				+":"+this.leftStart+","+this.leftEnd+","+this.rightStart+","+this.rightEnd;
		String out=this.chr+"\t"+this.start1+"\tINV"+i+"\tN\t<"+this.type+">\t"+qual+"\t"+filter+"\t"+info+"\t"+format+"\t"+formatResult+"\n";
		return out;
	}
	public String getchr(){
		return this.chr;
	}
	public int getLeftStart(){
		return this.leftStart;
	}
	public int getLeftEnd(){
		return this.leftEnd;
	}
	public int getRightStart(){
		return this.rightStart;
	}
	public int getRightEnd(){
		return this.rightEnd;
	}
	public boolean checkOrientation(){
		if(this.orientation[0]+this.orientation[1]>0 && this.orientation[2]+this.orientation[3]>0){
			return true;
		}else{
			return false;
		}
	}
	public void setReference(int[] refDepth){
		this.reference[0]=refDepth[0];
		this.reference[1]=refDepth[1];
		this.reference[2]=refDepth[2];
		this.reference[3]=refDepth[3];
		double[] gen=genotype(this.reference,this.orientation);
		//find major allele and then alt allele
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
		this.genotype=maxID;
		this.genotypeQuality=(int)(10*(max-mid));
	}
	public static double[] genotype(int[] ref,int[] variant){
		double[] out=new double[3];
		int refAll=ref[0]+ref[1]+ref[2]+ref[3];
		int varAll=variant[0]+variant[1]+variant[2]+variant[3];
		// A=10 B=5 15!/(10!*5!)=c is removed because every one has it.
		//update version 1.1 for genotype 0.25 0.5 0.25
		//P(G=AA)=P(G=BB)=0.25 //P(G=AB)=0.5
		out[0]=refAll*Math.log10(1-MergedSV.variantRate)+varAll*Math.log10(MergedSV.variantRate);
		out[0]+=Math.log10(0.25);//*0.25
		out[1]=(refAll+varAll)*Math.log10(0.5);
		out[1]+=Math.log10(0.5);
				//refAll*Math.log10((1-MergedSV.variantRate+MergedSV.variantRate)/2)+varAll*Math.log10((MergedSV.variantRate+1-MergedSV.variantRate)/2);
		out[2]=refAll*Math.log10(MergedSV.variantRate)+varAll*Math.log10(1-MergedSV.variantRate);
		out[2]+=Math.log10(0.25);
//		out[0]+=Math.log10(0.954529);
//		out[1]+=Math.log10(0.044942);
//		out[2]+=Math.log10(0.000529);
		return out;
	}
	public int[] getReference(){
		return this.reference;
	}
	public int[] getOrientation(){
		return this.orientation;
	}
	@Override
	public int compareTo(MergedSV o) {
		// TODO Auto-generated method stub
		int out=SV.compareChr(this.chr,o.chr);
		if(out==0){
			return this.start1-o.start1;
		}else{
			return out;
		}
	}
	public int getOrientationTotal(){
		return this.orientation[0]+this.orientation[1]+this.orientation[2]+this.orientation[3];
	}
	public static ArrayList<MergedSV> deletedRelate(ArrayList<MergedSV> mSV) {
		// TODO Auto-generated method stub
		boolean[] deleted=new boolean[mSV.size()];
		Arrays.fill(deleted, false);
		for(int i=0;i<mSV.size();i++){
			for(int j=i;j<mSV.size();j++){
				if(i!=j){
					if(overlap80(mSV.get(i),mSV.get(j))){
						deleted[i]=true;
					}
				}
			}
		}
		for(int i=mSV.size()-1;i>=0;i--){
			if(deleted[i]){
				mSV.remove(i);
			}
		}
		return mSV;
	}
	//if mSV1 is contained by mSV2 and mSV1>80%*mSV2 ; remove mSV1
	private static boolean overlap80(MergedSV mSV1, MergedSV mSV2) {
		// TODO Auto-generated method stub
		if(mSV1.start1<=mSV2.start2 && mSV2.start1<=mSV1.start2){
			double len1=mSV2.start2-mSV1.start1;
			double len2=mSV1.start2-mSV2.start1;
			if(mSV1.start2-mSV2.start1==0){
				return false;
			}
			double div=len1/len2;
			if( 0.666666<div && div<1.5){
				return true;
			}else{
				return false;
			}
		}else{
			return false;
		}
	}
}
