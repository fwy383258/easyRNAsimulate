package simulator;

import java.util.ArrayList;

public class Bed12 {
	private String chr=null;
	private int start=0;
	private int end=0;
	private String name=null;
	private double score=0.0;
	private char strand=' ';
	private int thick_start=0;
	private int thick_end=0;
	private int item_rgb=0;
	private int block_count=0;
	private ArrayList<Integer> block_sizes=null;
	private ArrayList<Integer> block_starts=null;
	private ArrayList<String> ids=null;
	
	public Bed12() {
		this.block_sizes = new ArrayList<>();
		this.block_starts = new ArrayList<>();
		this.ids = new ArrayList<>();
	}
	public Bed12(String chr, int start, int end, String name, double score, char strand, int thick_start,
			int thick_end, int item_rgb, int block_count, ArrayList<Integer> block_sizes, ArrayList<Integer> block_starts) {
		super();
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.name = name;
		this.score = score;
		this.strand = strand;
		this.thick_start = thick_start;
		this.thick_end = thick_end;
		this.item_rgb = item_rgb;
		this.block_count = block_count;
		this.block_sizes = block_sizes;
		this.block_starts = block_starts;
	}
	public String getChr() {
		return chr;
	}
	public void setChr(String chr) {
		this.chr = chr;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public char getStrand() {
		return strand;
	}
	public void setStrand(char strand) {
		this.strand = strand;
	}
	public int getThick_start() {
		return thick_start;
	}
	public void setThick_start(int thick_start) {
		this.thick_start = thick_start;
	}
	public int getThick_end() {
		return thick_end;
	}
	public void setThick_end(int thick_end) {
		this.thick_end = thick_end;
	}
	public int getItem_rgb() {
		return item_rgb;
	}
	public void setItem_rgb(int item_rgb) {
		this.item_rgb = item_rgb;
	}
	public int getBlock_count() {
		return block_count;
	}
	public void setBlock_count(int block_count) {
		this.block_count = block_count;
	}
	public ArrayList<Integer> getBlock_sizes() {
		return block_sizes;
	}
	public void setBlock_sizes(ArrayList<Integer> block_sizes) {
		this.block_sizes = block_sizes;
	}
	public ArrayList<Integer> getBlock_starts() {
		return block_starts;
	}
	public void setBlock_starts(ArrayList<Integer> block_starts) {
		this.block_starts = block_starts;
	}
	public ArrayList<String> getIds() {
		return ids;
	}
	public void setIds(ArrayList<String> ids) {
		this.ids = ids;
	}
	public static String getHeader() {
		return "chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
	}
    @Override
	public String toString() {
		StringBuffer out = new StringBuffer();
		out.append(chr);
		out.append('\t');
		out.append(start);
		out.append('\t');
		out.append(end);
		out.append('\t');
		out.append(name);
		out.append('\t');
		out.append(score);
		out.append('\t');
		out.append(strand);
		out.append('\t');
		out.append(thick_start);
		out.append('\t');
		out.append(thick_end);
		out.append('\t');
		out.append(item_rgb);
		out.append('\t');
		out.append(block_count);
		out.append('\t');
		out.append(blockSizesToString());
		out.append('\t');
		out.append(blockStartsToString());
		return out.toString();
	}
	public String blockSizesToString() {
		return aListToString(this.block_sizes, "0", ",");
	}
	public String blockStartsToString() {
		return aListToString(this.block_starts, "0", ",");
	}
	public String idsToName() {
		this.name = aListToString(this.ids, "Default", ",");
		return name;
	}
	public static <T> String aListToString(ArrayList<T> list, String null_out, String ifs) {
		if (list == null || list.size() < 1) {
			return null_out;
		}
		StringBuffer out = new StringBuffer();
		out.append(list.get(0));
		for (int i = 1; i < list.size(); i++) {
			out.append(ifs);
			out.append(list.get(i));
		}
		return out.toString();
	}
}
