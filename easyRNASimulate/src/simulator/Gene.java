package simulator;

import java.util.ArrayList;

public class Gene {

	private Chromosome chr=null;
	private String id=null;
	private ArrayList<Transcript> scripts=null;
	private int start=0;
	private int end=0;
	
	public Gene(){
		this.scripts = new ArrayList<>(); 
	}
	
	public Gene(Chromosome chr, String id, ArrayList<Transcript> scripts, int start, int end){
		super();
		this.chr = chr;
		this.id = id;
		this.scripts = scripts;
		this.start = start;
		this.end = end;
	}
	
	public Chromosome getChr() {
		return chr;
	}

	public void setChr(Chromosome chr) {
		this.chr = chr;
	}

	public String getId() {
		return id;
	}
	
	public void setId(String id) {
		this.id = id;
	}
	
	public ArrayList<Transcript> getScripts() {
		return scripts;
	}
	
	public Transcript getScript(int index) {
		return scripts.get(index);
	}
	
	public void addScript(Transcript script) {
		this.scripts.add(script);
	}
	
	public void setScripts(ArrayList<Transcript> scripts) {
		this.scripts = scripts;
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
	
}
