package inversion;

public class IRPair  implements Comparable<IRPair> {
	//public String chr;
	public int leftStart;
	public int leftEnd;
	public int rightStart;
	public int rightEnd;
	public IRPair(String line1,String line2){
		line1.replace("\n", "");
		String[] info=line1.split("\t| ");
		//this.chr=info[0];
		this.leftStart=Integer.parseInt(info[1]);
		this.leftEnd=Integer.parseInt(info[2]);
		line2.replace("\n", "");
		info=line2.split("\t| ");
		this.rightStart=Integer.parseInt(info[1]);
		this.rightEnd=Integer.parseInt(info[2]);
	}
	@Override
	public int compareTo(IRPair o) {
		//int out=SV.compareChr(this.chr,o.chr);
		int out=0;
		if(out==0){
			return this.leftStart-o.leftStart;
		}else{
			return out;
		}
	}
	public String toString(){
		return this.leftStart+"\t"+this.leftEnd+"\t"+this.rightStart+"\t"+this.rightEnd;
		//return this.chr+"\t"+this.leftStart+"\t"+this.leftEnd+"\t"+this.rightStart+"\t"+this.rightEnd;
	}
	public String toString(String chr){
		return chr+"\t"+this.leftStart+"\t"+this.leftEnd+"\t"+this.rightStart+"\t"+this.rightEnd;
		//return this.chr+"\t"+this.leftStart+"\t"+this.leftEnd+"\t"+this.rightStart+"\t"+this.rightEnd;
	}
}
