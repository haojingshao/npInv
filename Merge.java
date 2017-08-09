package inversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.sun.tools.javac.code.Attribute.Array;

public class Merge {
	public static final boolean test=Main.test;
	public static final int MergedSize=Main.mergedSize; // too large to include wrong alignment
	//don't write any global ArrayList
	//	private static ArrayList<MergedSV> finalSV=new ArrayList<MergedSV>();
	public static ArrayList<MergedSV> MergeSignal(ArrayList<SV> svList){
		HashMap<Integer, Boolean> analysized=new HashMap<Integer, Boolean>();
		ArrayList<MergedSV> finalSV=new ArrayList<MergedSV>();
		for(int i=0;i<svList.size();i++){
			if(!analysized.containsKey(i)){
				ArrayList<SV> svArray=new ArrayList<SV>();
				svArray.add(svList.get(i));
				if(test){System.out.print(i+"\t");}
				int j=i+1;// id for next value
				if(j<svList.size()){
					SV next=svList.get(j);
					while((!analysized.containsKey(j)) && (next.getStart1()-svList.get(i).getStart1()<=MergedSize) ){
						if(check(next,svList.get(i))){
							if(test){System.out.print(j+"\t"+!analysized.containsKey(j)+"\t");}
							svArray.add(svList.get(j));
							analysized.put(j, true);
							if(test){System.out.print(!analysized.containsKey(j)+"\t");}
						}
						j++;
						if (j<svList.size()){
							next=svList.get(j);
						}else{
							break;
						}
					}
				}
				MergedSV tmp=new MergedSV(svArray);
				if(test){System.out.print("\n");}
				finalSV.add(tmp);
			}else{
				//do nothing if contains
			}
		}
		return finalSV;
	}
	private static boolean check(SV next, SV sv) {
		return next.getchr().equals(sv.getchr()) && next.getChr2().equals(sv.getChr2()) &&
				  Math.abs(next.getStart2()-sv.getStart2())<MergedSize;
	}
	//for deletion Inversion
	//remove this function too much FDR
	public static ArrayList<MergedSV> MergeComplex(ArrayList<MergedSV> mSV) {
		//sort by end
		Collections.sort(mSV,new Comparator<MergedSV>() {
			@Override
			public int compare(MergedSV o1, MergedSV o2) {
				// TODO Auto-generated method stub
				int out=SV.compareChr(o1.chr,o2.chr);
				if(out==0){
					return o1.start2-o2.start2;
				}else{
					return out;
				}
			}	
		});
		//for start2
		int[] deletedID=new int[mSV.size()];
		Arrays.fill(deletedID,0);
		for(int i=0;i<mSV.size();i++){
			//select the read with the highest inversion reads
			for(int j=i+1;j<mSV.size() && mSV.get(j).start2-mSV.get(i).start2<MergedSize;j++){
				if(checkHalfSV(mSV.get(i),mSV.get(j))){
					mSV.add(new MergedSV(mSV.get(i),mSV.get(j)));//gain even more SV record!!!
					deletedID[i]++;
					deletedID[j]++;
				}
			}
		}
		for(int i=deletedID.length-1;i>=0;i--){
			//if(test){System.out.println("mSVsize= "+mSV.size()+" ;deletedID= "+i+";deletedSize"+deletedID.length);}
			if(deletedID[i]>0){			
				mSV.remove(i);
			}
		}
		//for start1
		Collections.sort(mSV);
		deletedID=new int[mSV.size()];
		Arrays.fill(deletedID,0);
		//if(test){System.out.println("2mSVsize= "+mSV.size()+" ;deletedSize= "+deletedID.size());}
		for(int i=0;i<mSV.size();i++){
			//select the read with the highest inversion reads
			for(int j=i+1;j<mSV.size() && mSV.get(j).start1-mSV.get(i).start1<MergedSize;j++){
				if(checkHalfSV(mSV.get(i),mSV.get(j))){
					mSV.add(new MergedSV(mSV.get(i),mSV.get(j)));
					deletedID[i]++;
					deletedID[j]++;
				}
			}
		}
		for(int i=deletedID.length-1;i>=0;i--){
			if(deletedID[i]>0){
				mSV.remove(i);
			}
		}
		Collections.sort(mSV);
		deletedID=new int[mSV.size()];
		Arrays.fill(deletedID,0);
		//if(test){System.out.println("3mSVsize= "+mSV.size()+" ;deletedSize= "+deletedID.size());}
		for(int i=0;i<mSV.size();i++){
			int j=1;
			while(i+j<mSV.size() && Math.abs(mSV.get(i).start1-mSV.get(i+j).start1)<MergedSize
					&& Math.abs(mSV.get(i).start2-mSV.get(i+j).start2)<MergedSize){
				
				if(mSV.get(i).getOrientationTotal()>mSV.get(i+j).getOrientationTotal()){
					deletedID[i+j]++;
					if(test){System.out.println("deleted sequence "+mSV.get(i+j).start1);}
					j++;
				}else{
					deletedID[i]++;
					if(test){System.out.println("deleted sequence "+mSV.get(i).start1);}
					break;
				}
			}
		}
		//if(test){System.out.println("4mSVsize= "+mSV.size()+" ;deletedSize= "+deletedID.size());}
		for(int i=deletedID.length-1;i>=0;i--){
			if(deletedID[i]>0){
				mSV.remove(i);
			}
		}
		//if(test){System.out.println("5mSVsize= "+mSV.size()+" ;deletedSize= "+deletedID.size());}
		return mSV;		
	}
	private static boolean checkHalfSV(MergedSV MergedSV, MergedSV MergedSV2) {
		if((MergedSV.orientation[0]+MergedSV.orientation[1]>0 && MergedSV.orientation[2]+MergedSV.orientation[3]==0
		&& MergedSV2.orientation[2]+MergedSV2.orientation[3]>0 && MergedSV2.orientation[0]+MergedSV2.orientation[1]==0)||
		(MergedSV.orientation[2]+MergedSV.orientation[3]>0 && MergedSV.orientation[0]+MergedSV.orientation[1]==0
		&& MergedSV2.orientation[0]+MergedSV2.orientation[1]>0 && MergedSV2.orientation[2]+MergedSV2.orientation[3]==0)
				){
			return true;
		}else{
			return false;
		}
	}
	//this is new function
	//if there is an IR database ,1 group all sv signal into IR 2 then summary the remaining sv
	//else run 2 summary the sv
	public static ArrayList<MergedSV> MergeSignal(ArrayList<SV> sv, String irFile, String chr) {
		if(irFile.equals("none") ){
			return MergeSignal(sv);
		}else {
			File tmpFile = new File(irFile);
			if(tmpFile.exists() && !tmpFile.isDirectory()) { 
				ArrayList<IRPair> IR=readirFile(irFile,chr); 
				ReturnMergeIRSignal returnMergeIRSignal=null;
				returnMergeIRSignal=MergeIRSignal(sv,IR,chr);// first return to finalSV
				returnMergeIRSignal.finalSV.addAll(MergeSignal(returnMergeIRSignal.sv));// then add to finalSV
				Collections.sort(returnMergeIRSignal.finalSV); //sort two result
				return returnMergeIRSignal.finalSV;
			}else{				
			    System.out.println("IR file not exists. Running without it\n");
			    return MergeSignal(sv);
			}
		}
	}
	private static ReturnMergeIRSignal MergeIRSignal(ArrayList<SV> sv, ArrayList<IRPair> IR, String chr2) {
		//split chromosome and O(n2)
		//ArrayList<String> chrAll=new ArrayList<String>();
		//Map<chr,<IRfirst,IRend,SVfirst,SVend>>
		//They don't need to split chr by yourself! they don't have chromosome
		//Map<String,Integer[]> splitChr=getChrStartEnd(sv,IR);
		//for each chr
		int[] deletedID=new int[sv.size()];
		Arrays.fill(deletedID, 0);
		ArrayList<MergedSV> finalSV =new ArrayList<MergedSV>() ;
		for (int iIR=0;iIR<IR.size();iIR++){
			ArrayList<SV> group=new ArrayList<SV>();
			for(int iSV=0;iSV<sv.size();iSV++){
				if(checkIRSV(IR.get(iIR),sv.get(iSV))){
					group.add(sv.get(iSV));
					deletedID[iSV]++;			
				}
			}
			//System.out.println("groupSize="+group.size()+" IR= "+IR.get(iIR).toString());
			if(group.size()>0){
				finalSV.add(new MergedSV(group,IR.get(iIR)));
			}
		}
		for (int i=deletedID.length-1;i>=0;i--){
			if(deletedID[i]>0){
				sv.remove(i);
			}
		}
		ReturnMergeIRSignal returnMergeIRSignal=new ReturnMergeIRSignal(sv,finalSV);
		return returnMergeIRSignal;
	}
	//deleted function
	/*
	public static Map<String,Integer[]> getChrStartEnd(ArrayList<SV> sv, ArrayList<IRPair> IR){
		Map<String,Integer[]> splitChr=new HashMap<String,Integer[]>();
		String chr="";
		String thisChr="";
		Integer[] tmpInt=new Integer[4];
		Arrays.fill(tmpInt, 0);
		//I should somewhat reuse this code in a function??
		for(int i=0;i<IR.size();i++){
			thisChr=IR.get(i).chr;
			if(!thisChr.equals(chr)){
				tmpInt=new Integer[4];
				Arrays.fill(tmpInt, 0);
				tmpInt[0]=i;
				tmpInt[1]=i;// in case there is nothing
				splitChr.put(thisChr, tmpInt);
				chr=thisChr;
			}else{
				splitChr.get(thisChr)[1]=i;
			}
		}
		chr="";
		for(int i=0;i<sv.size();i++){
			thisChr=sv.get(i).chr;
			if(!splitChr.containsKey(thisChr)){
				tmpInt=new Integer[4];
				Arrays.fill(tmpInt, 0);
				splitChr.put(thisChr, tmpInt);
			}
			if(!thisChr.equals(chr)){
				splitChr.get(thisChr)[2]=i;
				splitChr.get(thisChr)[3]=i;
				chr=thisChr;
			}else{
				splitChr.get(thisChr)[3]=i;
			}
		}
		return splitChr;
	}
	*/
	private static boolean checkIRSV(IRPair ir, SV sv) {
		if(ir.leftStart-MergedSize<sv.start1 && sv.start1<ir.leftEnd+MergedSize
				&& ir.rightStart-MergedSize<sv.start2 && sv.start2<ir.rightEnd+MergedSize){
			return true;
		}else{
			return false;
		}
	}
	private static ArrayList<IRPair> readirFile(String irFile,String chr) {
		// TODO Auto-generated method stub
		ArrayList<IRPair> IR=new ArrayList<IRPair>();
		try (BufferedReader br = new BufferedReader(new FileReader(irFile))) {
		    String line1;
		    String line2;
		    while ((line1 = br.readLine()) != null) {
		    		String [] split=line1.split("\t| ");// \t and space
		    		if(split[0].equals(chr)){
			    		line2=br.readLine();
			    		IR.add(new IRPair(line1,line2));
		    		}else{
		    			//don't add 
		    			line2=br.readLine();
		    		}
		    }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Collections.sort(IR);
		return IR;
	}


}
