package simulator;

import java.util.ArrayList;
import java.util.HashSet;

public class Method {
	
	static char[] bases = {'A', 'G', 'T', 'C'};
	
	static <T> void getNoneRepeat(ArrayList<T> list, int count, HashSet<T> out) {
		if (list.size() > count) {
			while(out.size() < count) {
				out.add(list.get(Method.randIntReach(0, list.size() - 1)));
			}
		}
		else {
			System.out.println("Warning: Not enough elements for picking");
		}
	}
	
	static <T> T getRandElement(ArrayList<T> list) {
		return list.get(Method.randIntReach(0, list.size() - 1));
	}
	
	static int randIntReach(int low, int high) {
		int out = 0;
		if(high >= low) {
			out = (int) (Math.random() * (high - low + 1)) + low;
		}
		else {
			System.out.println("Warning: " + low + " is greater than " + high);
		}
		return out;
	}
	
	static int[] divideInto(int total, int size) {
		if (size <= 0) {
			return null;
		}
		else if (size == 1) {
			int[] out = {total};
			return out;
		}
		else {
			int[] out = new int[size];
			ArrayList<Integer> temp = new ArrayList<Integer>();
			for (int i = 1; i < size; i++) {
				temp.add((int) (Math.random() * (total + 1)));
			}
			temp.sort(null);
			out[0] = temp.get(0);
			for (int i = 1; i < size - 1; i++) {
				out[i] = temp.get(i) - temp.get(i - 1);
			}
			out[size - 1] = total - temp.get(size - 2);
			return out;
		}
	}
	
	/*
	 * search the first integer not less than target number in list
	 * return the index of it
	 */
	static int searchMinNoLess(int target, ArrayList<Integer> inc_seq) {
		int out = -1;
		int l = 0;
		int r = inc_seq.size() - 1;
		if (target <= inc_seq.get(0)) {
			out = 0;
		}
		else if(target <= inc_seq.get(r)){
			int m = 0;
			while (l < r) {
				m = (l + r) >> 1;
				if (l == m) {
					out = r;
					break;
				}
				if (target < inc_seq.get(m)) {
					r = m;
				}
				else if (target > inc_seq.get(m)) {
					l = m;
				}
				else {
					while (inc_seq.get(m) == target) {
						out = m;
						m--;
					}
					break;
				}
			}
		}
		return out;
	}
	
	static void randomErrorSeq(StringBuffer seq) {
		for (int i = 0 ; i < seq.length(); i++) {
			double rand = Math.random();
			if (rand < 0.01) {
				if (rand < 0.0002) {
					if (rand < 0.0001) {
						seq.deleteCharAt(i);
					}
					else {
						seq.insert(i, randBase('N'));
					}
				}
				else {
					seq.setCharAt(i, randBase(seq.charAt(i)));
				}
			}
		}
	}
	
	static char randBase(char old) {
		char out = old;
		while(out == old) {
			out = bases[randIntReach(0, bases.length - 1)];
		}
		return out;
	}
}
