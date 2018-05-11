package simulator;

import java.util.ArrayList;
import java.util.Collections;

public class Transcript {

	private Gene gene=null;
	private String id=null;
	private ArrayList<Exon> exons=null;
	private int start=0;
	private int end=0;
	private int base_sum=0;
	private char strand=' ';
	private boolean circ_flag=false;
	
	public Transcript(){
		this.exons = new ArrayList<>();
	}
	
	public Transcript(Gene gene, String id, ArrayList<Exon> exons, int start, int end, int base_sum, char strand, boolean circ_flag) {
		super();
		this.gene = gene;
		this.id = id;
		this.exons = exons;
		this.start = start;
		this.end = end;
		this.base_sum = base_sum;
		this.strand = strand;
		this.circ_flag = circ_flag;
	}

	public Gene getGene() {
		return gene;
	}
	
	public void setGene(Gene gene) {
		this.gene = gene;
	}
	
	public String getId() {
		return id;
	}
	
	public void setId(String id) {
		this.id = id;
	}
	
	public ArrayList<Exon> getExons() {
		return exons;
	}
	
	public Exon getExon(int index) {
		return exons.get(index);
	}
	
	public void addExon(Exon exon) {
		this.exons.add(exon);
	}
	
	public void setExons(ArrayList<Exon> exons) {
		this.exons = exons;
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

	public int getBase_sum() {
		return base_sum;
	}

	public void setBase_sum(int base_sum) {
		this.base_sum = base_sum;
	}

	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}

	public boolean isCirc_flag() {
		return circ_flag;
	}

	public void setCirc_flag(boolean circ_flag) {
		this.circ_flag = circ_flag;
	}

	/*
	 * param:less_front means sort from less to greater while true
	 */
	public void sortExons(boolean less_front) {
		this.quickSortExons(0, this.exons.size()- 1);
		if (!less_front) {
			Collections.reverse(this.exons);
		}
	}
	
	private void quickSortExons(int start, int end) {
		if (end - start <= 8) {
			this.insertSortExons(start, end);
			return;
		}
		int left = start;
		int right = end;
		int middle = (left + right) >> 1;
		int key = this.exons.get(middle).getStart();
		
		while(left < right) {
			while (this.exons.get(left).getStart() <= key) {
				left++;
			}
			while (this.exons.get(right).getStart() >= key) {
				right--;
			}
			if (left < right) {
				this.exons.set(right, this.exons.set(left, this.exons.get(right)));
			}
			else if (left < (start+end) >> 1) {
				this.exons.set(left, this.exons.set(middle, this.exons.get(left)));
				right = left;
			}
			else if (right > middle){
				this.exons.set(right, this.exons.set(middle, this.exons.get(right)));
				left = right;
			}
		}
		
		this.quickSortExons(start, left - 1);
		this.quickSortExons(right + 1, end);
	}
	
	private void insertSortExons(int start, int end) {
		for (int i=start + 1; i <= end; i++) {
			int key = this.exons.get(i).getStart();
			for (int j=start; j < i; j++) {
				if (key < this.exons.get(j).getStart()) {
					this.exons.add(j, this.exons.get(i));
					this.exons.remove(i + 1);
					break;
				}
			}
		}
	}
	
}
