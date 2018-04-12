package simulator;

import java.util.ArrayList;

public class Exon {
	
	private Transcript script=null;
	private String base_seq=null;
	private String chr_symbol=null;
	private int chr_num=0;
	private int start=0;
	private int end=0;
	
	
	public Transcript getScript() {
		return script;
	}

	public void setScript(Transcript script) {
		this.script = script;
	}

	String getBase_seq() {
		return base_seq;
	}

	void setBase_seq(String base_seq) {
		this.base_seq = base_seq;
	}

	String getChr_symbol() {
		return chr_symbol;
	}
	
	void setChr_symbol(String chr_symbol) {
		this.chr_symbol = chr_symbol;
		this.chr_num = this.chrSymbolToNum(chr_symbol);
	}
	
	int getChr_num() {
		return chr_num;
	}
	
	void setChr_num(int chr_num) {
		this.chr_num = chr_num;
		this.chr_symbol = this.chrNumToSymbol(chr_num);
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
	
	public int chrSymbolToNum(String chr_symbol) {
		int out = -1;
		if (chr_symbol.length() > 3 && chr_symbol.length() < 6) {
			String temp = chr_symbol.substring(3);
			if (temp.equals("X")) {
				out = 23;
			}
			else if (temp.equals("Y")) {
				out = 24;
			}
			else if (temp.equals("M")) {
				out = 25;
			}
			else {
				out = Integer.parseInt(temp);
			}
		}
		return out;
	}
	
	public String chrNumToSymbol(int chr_num) {
		String out = null;
		if (chr_num > 25 || chr_num <= 0) {
			 System.out.println("Undefine CHR");
		}
		else if (chr_num==23) {
			out = "chrX";
		}
		else if (chr_num==24) {
			out = "chrY";
		}
		else if (chr_num==25) {
			out = "chrM";
		}
		else {
			out = "chr" + chr_num;
		}
		return out;
	}
	
	int searchMinNoLess(int target, ArrayList<Exon> inc_seq) {
		int out = -1;
		int l = 0;
		int r = inc_seq.size() - 1;
		if (inc_seq.size() == 0||target <= inc_seq.get(0).getStart()) {
			out = 0;
		}
		else if(target <= inc_seq.get(r).getStart()){
			int m = 0;
			while (l < r) {
				m = (l + r) >> 1;
				if (l == m) {
					out = r;
					break;
				}
				if (target < inc_seq.get(m).getStart()) {
					r = m;
				}
				else if (target > inc_seq.get(m).getStart()) {
					l = m;
				}
				else {
					while (target == inc_seq.get(m).getStart()) {
						out = m;
						m--;
					}
					break;
				}
			}
		}
		else {
			out = r + 1;
		}
		return out;
	}
	
}
