package simulator;

import java.util.ArrayList;
import java.util.HashMap;

public class Chromosome {
	
	private String id=null;
	private ArrayList<Gene> genes=null;
	private ArrayList<Integer> all_junctions=null;
	private HashMap<Integer, ArrayList<BackJunction>> junction_map=null;
	
	public Chromosome(){
		this.genes = new ArrayList<>();
		this.all_junctions = new ArrayList<>();
		this.junction_map = new HashMap<>();
	}
	
	public Chromosome(String id, ArrayList<Gene> genes, ArrayList<Integer> all_junctions, HashMap<Integer, ArrayList<BackJunction>> junction_map){
		super();
		this.id = id;
		this.genes = genes;
		this.all_junctions = all_junctions;
		this.junction_map = junction_map;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public ArrayList<Gene> getGenes() {
		return genes;
	}

	public Gene getGene(int index) {
		return genes.get(index);
	}
	
	public void addGene(Gene gene) {
		this.genes.add(gene);
	}
	
	public void setGenes(ArrayList<Gene> genes) {
		this.genes = genes;
	}

	public ArrayList<Integer> getAll_junctions() {
		return all_junctions;
	}

	public void setAll_junctions(ArrayList<Integer> all_junctions) {
		this.all_junctions = all_junctions;
	}

	public HashMap<Integer, ArrayList<BackJunction>> getJunction_map() {
		return junction_map;
	}

	public void setJunction_map(HashMap<Integer, ArrayList<BackJunction>> junction_map) {
		this.junction_map = junction_map;
	}

	public static int chrSymbolToNum(String chr_symbol) {
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
	
	public static String chrNumToSymbol(int chr_num) {
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
	
}
