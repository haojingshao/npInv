package inversion;

public class startEnd {
	public int start;
	public int end;
	public boolean negativeStrandFlag;
	public startEnd(int alignmentStart, int alignmentEnd,boolean strand) {
		// TODO Auto-generated constructor stub
		this.start=alignmentStart;
		this.end=alignmentEnd;
		this.negativeStrandFlag=strand;
	}
	@Override
	public String toString(){
		//override
		return this.start+" "+this.end+" "+this.negativeStrandFlag;
	}
}
