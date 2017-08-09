package inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


// read depth for every SV
public class RefDepth {
	// unused function... deleted
	public static ArrayList<MergedSV> getDepth(ArrayList<MergedSV> mSV, String input) throws IOException {
		// TODO Auto-generated method stubz
		for(int i=0;i<mSV.size();i++){
			//ls= left start
			int ls=mSV.get(i).getLeftStart();
			int le=mSV.get(i).getLeftEnd();
			int rs=mSV.get(i).getRightStart();
			int re=mSV.get(i).getRightEnd();
			//0= lf 1=lr 2=rf 3=rr
			int[] RefDepth=new int[4];
			RefDepth[0]=0;
			RefDepth[1]=0;
			RefDepth[2]=0;
			RefDepth[3]=0;
			//re-read a file many times
			final SamReaderFactory factory =
			         SamReaderFactory.makeDefault()
		             .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
		             SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		             .validationStringency(ValidationStringency.SILENT);

			final SamReader samReader = factory.open(new File(input));
			SAMRecord sam;
			SAMRecordIterator iter;
			//for left breakpoint
			iter=samReader.query(mSV.get(i).getchr(), ls ,le , true);
			while(iter.hasNext()){
				sam=iter.next();
				int start=sam.getAlignmentStart();
				int end=sam.getAlignmentEnd();
				if(sam.getReadNegativeStrandFlag()){
					//-
					if(start<ls && le<end){
						RefDepth[1]+=1;
					}
				}else{
					//+
					if(start<ls && le<end){
						RefDepth[0]+=1;
					}
				}
			}
			//for right break point			
			iter=samReader.query(mSV.get(i).getchr(), rs ,re , true);
			while(iter.hasNext()){
				sam=iter.next();
				int start=sam.getAlignmentStart();
				int end=sam.getAlignmentEnd();
				if(sam.getReadNegativeStrandFlag()){
					//-
					if(start<rs && re<end){
						RefDepth[3]+=1;
					}
				}else{
					//+
					if(start<rs && re<end){
						RefDepth[2]+=1;
					}
				}
			}
			mSV.get(i).setReference(RefDepth);
			samReader.close();
		}
		
		return mSV;
	}
	public static ArrayList<MergedSV> removeLowCount(ArrayList<MergedSV> mSV) {
		Iterator<MergedSV> iter=mSV.iterator();
		while(iter.hasNext()){
			MergedSV tmp=iter.next();
			int[] inv=tmp.getOrientation();
			if(inv[0]+inv[1]>1 && inv[2]+inv[3]>1){
			}else{
				iter.remove();
			}
		}
		return mSV;
	}

}
