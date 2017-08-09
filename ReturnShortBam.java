package inversion;

public class ReturnShortBam {
	public startEnd[] depth;
	public ClipRead[] signal;
	public String chr;
	public ReturnShortBam(startEnd[] depthShort, ClipRead[] clip, String chr2) {
		this.depth=depthShort;
		this.signal=clip;
		this.chr=chr2;
	}
}
