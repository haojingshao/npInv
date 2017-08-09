package inversion;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class ClipRead {
	public static int mergedClipRead=Main.mergedClipRead;
	public static int minimumClipRead=Main.minimumClipRead;
	public int pos;
	public boolean forwardOrNot;
	public int count;
	//for single
	public ClipRead(int alignmentStart, boolean b) {
		// TODO Auto-generated constructor stub
	//	this.chr=referenceName;
	//	this.chr=0;
		this.pos=alignmentStart;
		this.forwardOrNot=b;
		this.count=1;
	}
	//for after merge
	public ClipRead(ArrayList<ClipRead> group) {
//		this.chr=group.get(0).chr;
	//	this.chr=0;
		this.forwardOrNot=group.get(0).forwardOrNot;
		this.pos=group.get(group.size()/2).pos;//median pos
		this.count=group.size();
	}
	public static ClipRead[] mergeAll(ClipRead[] clipRaw) {
		// TODO Auto-generated method stub
		ArrayList<ClipRead> result=new ArrayList<ClipRead>();
		boolean[] notAnalyzed=new boolean[clipRaw.length];
		Arrays.fill(notAnalyzed, true);
		for(int i=0;i<clipRaw.length;i++){
			notAnalyzed[i]=false;
			int offset=i+1;
			ArrayList<ClipRead> group=new ArrayList<ClipRead>();
			group.add(clipRaw[i]);
			while(offset<clipRaw.length && notAnalyzed[offset]){
				//merge to the first one OR merge relat
				//version 0.0 have bug with close inversion
				if(clipRaw[offset].pos-clipRaw[i].pos<mergedClipRead
					&& !clipRaw[offset].forwardOrNot^clipRaw[i].forwardOrNot){
					group.add(clipRaw[offset]);
					notAnalyzed[offset]=false;
				}
				offset++;
			}
			if(group.size()>=minimumClipRead){
				result.add(new ClipRead(group));
			}
		}
		ClipRead[] resultArray=new ClipRead[result.size()];
		resultArray=result.toArray(resultArray);
		return resultArray;
	}
	public String toString(){
		return this.pos+"\t"+this.forwardOrNot+"\t"+this.count;
	}
	public static ArrayList<MergedSV> mSVwithClipRead(ArrayList<MergedSV> mSV, ClipRead[] clip, startEnd[] depth) {
		// TODO Auto-generated method stub
		ArrayList<MergedSV> newMSV=new ArrayList<MergedSV>();
		for(int iSV=0;iSV<mSV.size();iSV++){
			///only check the mSV with less homology size Or short SV
			MergedSV tmpMSV=mSV.get(iSV);
			if(tmpMSV.homo<=200 || tmpMSV.start2-tmpMSV.start1<=3000){
				MergedSV clipMSV=checkRegion(tmpMSV,clip);
				if(clipMSV.KeepOrNot){
					Main.getDepthSingle(clipMSV, depth);// already modified
					newMSV.add(clipMSV);
				}
			}
			
		}
		return newMSV;
	}
	private static MergedSV checkRegion(MergedSV tmpMSV, ClipRead[] clip) {
		// TODO Auto-generated method stub
		ArrayList<ClipRead> tmpClip1=new ArrayList<ClipRead>();
		ArrayList<ClipRead> tmpClip2=new ArrayList<ClipRead>();
		for(int i=0;i<clip.length;i++){
			if(tmpMSV.leftStart<=clip[i].pos && clip[i].pos<=tmpMSV.leftEnd){
				tmpClip1.add(clip[i]);
			}
			if(tmpMSV.rightStart<=clip[i].pos && clip[i].pos<=tmpMSV.rightEnd){
				tmpClip2.add(clip[i]);
			}
		}
		return new MergedSV(tmpMSV,tmpClip1,tmpClip2);
	}
	//deleted function
	/*
	private static Map<String, Integer[]> getChrStartEnd(ClipRead[] clip) {
		// TODO Auto-generated method stub
		Map<String,Integer[]> splitChr=new HashMap<String,Integer[]>();
		String chr="";
		String thisChr="";
		//int chr=0;
		//int thisChr=0;
		Integer[] tmpInt=new Integer[2];
		Arrays.fill(tmpInt, 0);
		//I should somewhat reuse this code in a function??
		for(int i=0;i<clip.length;i++){
			thisChr=clip[i].chr;
			if(!thisChr.equals(chr)){
			//if(thisChr == chr){
				tmpInt=new Integer[2];
				Arrays.fill(tmpInt, 0);
				tmpInt[0]=i;
				tmpInt[1]=i;// in case there is nothing
				splitChr.put(thisChr, tmpInt);
				chr=thisChr;
			}else{
				splitChr.get(thisChr)[1]=i;
			}
		}
		return splitChr;
	}
	*/
}
