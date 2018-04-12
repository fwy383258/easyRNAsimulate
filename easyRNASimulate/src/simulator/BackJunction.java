package simulator;

import java.util.ArrayList;

public class BackJunction {
	
	private Chromosome chr=null;
	private int start=0;
	private int end=0;
	private String start_seq = null;
	private String end_seq =null;
	private ArrayList<Transcript> start_scripts=null;
	private ArrayList<Transcript> end_scripts=null;
	private ArrayList<Transcript> both_scripts=null;

	
	public BackJunction() {
		this.start_scripts = new ArrayList<>();
		this.end_scripts = new ArrayList<>();
		this.both_scripts = new ArrayList<>();
	}

	public Chromosome getChr() {
		return chr;
	}

	public void setChr(Chromosome chr) {
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

	public String getStart_seq() {
		return start_seq;
	}

	public void setStart_seq(String start_seq) {
		this.start_seq = start_seq;
	}

	public String getEnd_seq() {
		return end_seq;
	}

	public void setEnd_seq(String end_seq) {
		this.end_seq = end_seq;
	}

	public ArrayList<Transcript> getStart_scripts() {
		return start_scripts;
	}
	
	public void addStart_script(Transcript script) {
		this.start_scripts.add(script);
	}

	public void setStart_scripts(ArrayList<Transcript> start_scripts) {
		this.start_scripts = start_scripts;
	}

	public ArrayList<Transcript> getEnd_scripts() {
		return end_scripts;
	}

	public void addEnd_script(Transcript script) {
		this.end_scripts.add(script);
	}
	
	public void setEnd_scripts(ArrayList<Transcript> end_scripts) {
		this.end_scripts = end_scripts;
	}
	public ArrayList<Transcript> getBoth_scripts() {
		return both_scripts;
	}

	public void addBoth_script(Transcript script) {
		this.both_scripts.add(script);
	}
	
	public void setBoth_scripts(ArrayList<Transcript> both_scripts) {
		this.both_scripts = both_scripts;
	}

}
