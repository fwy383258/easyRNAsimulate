package simulator;

import java.util.ArrayList;

public class Read {
	private String id=null;
	private ArrayList<Integer> positions=null;
	private boolean circ_flag=false;
	private boolean peak_flag=false;
	
	public Read() {
		this.positions = new ArrayList<>();
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public ArrayList<Integer> getPositions() {
		return positions;
	}

	public void setPositions(ArrayList<Integer> positions) {
		this.positions = positions;
	}

	public boolean isCirc_flag() {
		return circ_flag;
	}

	public void setCirc_flag(boolean circ_flag) {
		this.circ_flag = circ_flag;
	}

	public boolean isPeak_flag() {
		return peak_flag;
	}

	public void setPeak_flag(boolean peak_flag) {
		this.peak_flag = peak_flag;
	}
}
