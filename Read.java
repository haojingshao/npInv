package inversion;

import java.util.ArrayList;
import java.util.Collections;

public class Read {
//	private int mappingStart;
//	private int mappingEnd;
	private ArrayList<Alignment> AlignmentArray;
	public Read(ArrayList<Alignment> Alignment,int minSize){
		this.AlignmentArray=new ArrayList<Alignment>(minimumOverlap(Alignment,minSize));
//		if(true){
//			for(int i=0;i<this.AlignmentArray.size();i++){
//		//		System.out.println("s="+this.AlignmentArray.get(i).getMappingStart()+"\te="+this.AlignmentArray.get(i).getMappingEnd() );
//			}
//		}
//		int start=0;
//		int end=0;
//		for(int i=0;i<this.AlignmentArray.size();i++){
//			if(start==0 || this.AlignmentArray.get(i).getMappingStart()<start){
//				start=this.AlignmentArray.get(i).getMappingStart();
//			}
//			if(end==0 || end<this.AlignmentArray.get(i).getMappingEnd()){
//				end=this.AlignmentArray.get(i).getMappingEnd();
//			}
//		}
//		this.mappingStart=start;
//		this.mappingEnd=end;
	}
	public ArrayList<Alignment> getAlignmentArray(){
		return this.AlignmentArray;
	}
	public static ArrayList<Alignment> minimumOverlap(ArrayList<Alignment> startEnd,int minSize) {
		// time O(n**2)
		ArrayList<Integer> remainAlignment = new ArrayList<Integer>();
		boolean keep;
		for (int i=0;i<startEnd.size();i++){
			keep=true;
			if(startEnd.get(i).getMappingEnd()-startEnd.get(i).getMappingStart()<minSize){
				keep=false;//too short remove record
			}else{
				for (int j=0;j<startEnd.size();j++){
					if(j!=i){
						if((startEnd.get(j).getMappingStart()<=startEnd.get(i).getMappingStart() && startEnd.get(i).getMappingEnd()<startEnd.get(j).getMappingEnd() )||
								(startEnd.get(j).getMappingStart()<startEnd.get(i).getMappingStart() && startEnd.get(i).getMappingEnd()<=startEnd.get(j).getMappingEnd())){ // double largest ???
							keep=false;
							break;
						}
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
				if(startEnd.get(remainAlignment.get(i)).getMappingStart()==startEnd.get(remainAlignment.get(j)).getMappingStart() 
						&& startEnd.get(remainAlignment.get(i)).getMappingEnd()==startEnd.get(remainAlignment.get(j)).getMappingEnd()){
					remainAlignment.remove(i);
				}
			}
		}
		// TODO Auto-generated method stub
		ArrayList<Alignment> out=new ArrayList<Alignment>();
		for(int i=0;i<remainAlignment.size();i++){
			out.add(startEnd.get(remainAlignment.get(i)));
		}
		Collections.sort(out);
		return out;
	}
	//@override
}
