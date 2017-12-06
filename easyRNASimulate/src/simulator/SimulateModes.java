package simulator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

public class SimulateModes {

	private static int total_modes=1;
	
	static void run(ReadInputs args) {
		ArrayList<String> out_reads = null;
		SimulateModes sim = new SimulateModes();
		int mode = args.getSim_mode();
		switch (mode) {
		case 1:
			out_reads = sim.mode1Uniform(args);
			sim.writeFasta(args.getOut_file(), out_reads);
			break;
		case 2:
			break;
		default:
			break;
		}
		
		
	}
	
	ArrayList<String> mode1Uniform(ReadInputs args){
		ArrayList<String> out = new ArrayList<String>();
		ArrayList<ArrayList<ExonInfo>> exons = this.loadFiles(args.getRef_exon_gtf(), args.getRef_genome_file());
		
		int read_count = args.getRead_count();
		int map_rand = (int) (read_count * (1.0 - args.getSplice_input_scale()));
		int circ_rand = read_count - map_rand;
		int splice_rand = (int) (circ_rand * (1.0 - args.getCirc_scale()));
		circ_rand -= splice_rand;
		int[] chr_rand = this.divideInto(map_rand, exons.size());
		int[] chr_srand = this.divideInto(splice_rand, exons.size());
		int[] chr_crand = this.divideInto(splice_rand, exons.size());
		for (int i = 0; i < exons.size(); i++) {
			ArrayList<Integer> exon_seq = new ArrayList<Integer>();
			int random_max = this.getExonRandomMax(exons.get(i), args.getRead_length() - 1, exon_seq);
			while (chr_rand[i] > 0) {
				int offset = (int) (random_max * Math.random());
				int exon_num = this.searchMinNoLess(offset + 1, exon_seq);
				if (exon_num > 0) {
//					System.out.printf("%d %d %d ", exon_seq.get(exon_num - 1), exon_seq.get(exon_num),offset);
					offset -= exon_seq.get(exon_num - 1);
				}
//				System.out.println(offset);
				String read = exons.get(i).get(exon_num).getBase_seq().substring(offset, offset + args.getRead_length());
				StringBuffer temp_read = new StringBuffer();
				temp_read.append('>');
				temp_read.append(exons.get(i).get(exon_num).getChr_symbol());
				temp_read.append(':');
				temp_read.append(exons.get(i).get(exon_num).getStart_position() + offset);
				out.add(temp_read.toString());
				out.add(read);
				chr_rand[i] --;
			}
			exon_seq = new ArrayList<Integer>();
			random_max = this.getExonForwardSpliceRandomMax(exons.get(i), args.getMax_splice(), args.getMin_splice(), exon_seq);
			while (chr_srand[i] > 0) {
				int offset = (int) (random_max * Math.random()) + 1;
				int exon_1 = this.searchMinNoLess(offset, exon_seq);
				offset += exons.get(i).get(exon_1).getBase_seq().length() - exon_seq.get(exon_1) - args.getMin_splice();
				StringBuffer temp = new StringBuffer();
				temp.append('>');
				String read = this.getBoundSplice(exons.get(i), exon_1, offset, args.getRead_length(), temp);
//				String read = this.getBoundSplicePieces(exons.get(i), exon_1, offset, args.getRead_length(), temp, out);
				out.add(temp.toString());
				out.add(read);
				chr_srand[i] --;
			}
			while (chr_crand[i] > 0) {
				int offset = (int) (random_max * Math.random()) + 1;
				int exon_1 = this.searchMinNoLess(offset, exon_seq);
				offset += exons.get(i).get(exon_1).getBase_seq().length() - exon_seq.get(exon_1) - args.getMin_splice();
				StringBuffer temp = new StringBuffer();
				temp.append('>');
				String read = this.getBoundCirc(exons.get(i), exon_1, offset, args.getRead_length(), temp);
//				String read = this.getBoundSplicePieces(exons.get(i), exon_1, offset, args.getRead_length(), temp, out);
				out.add(temp.toString());
				out.add(read);
				chr_crand[i] --;
			}
		}
		return out;
	}
	
	ArrayList<ArrayList<ExonInfo>> loadFiles(String exon_gtf_file, String ref_genome_file){
		ArrayList<ArrayList<ExonInfo>> out = new ArrayList<ArrayList<ExonInfo>>();
		for (int i = 0; i < 25; i++) {
			ArrayList<ExonInfo> temp = new ArrayList<ExonInfo>();
			out.add(temp);
		}
		FileRead file_read = new FileRead();
		String temp_line = null;
		ExonInfo exon = null;

		
		BufferedReader exon_gtf = file_read.openFile(exon_gtf_file);
		try {
			while((temp_line = exon_gtf.readLine()) != null) {
				String[] cols = temp_line.split("\t");
				if (cols[2].equals("exon")) {
					int exon_start = Integer.parseInt(cols[3]);
					int exon_end = Integer.parseInt(cols[4]);
					exon = new ExonInfo();
					exon.setChr_symbol(cols[0]);
					exon.setStart_position(exon_start);
					exon.setEnd_position(exon_end);
					int chr_array = exon.getChr_num() - 1;
					if (chr_array >= 0 && chr_array < 25) {
						ArrayList<ExonInfo> exon_list = out.get(chr_array);
						int index = exon.searchMinNoLess(exon_start, exon_list);
						if (exon_list.size() > index && exon_start == exon_list.get(index).getStart_position() && exon_end == exon_list.get(index).getEnd_position()) {
							continue;
						}
						exon_list.add(index, exon);
					}
				}
			}
			file_read.closeFile(exon_gtf);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		BufferedReader genome_fasta = file_read.openFile(ref_genome_file);
		try {
			int chr_array = -1;
			while((temp_line = genome_fasta.readLine()) != null) {
				if (temp_line.charAt(0) == '>') {
					chr_array = exon.chrSymbolToNum(temp_line.substring(1)) - 1;
				}
				else if (chr_array >= 0) {
					ArrayList<ExonInfo> exon_list = out.get(chr_array);
					for (int i = 0; i < exon_list.size(); i++) {
						exon = exon_list.get(i);
						exon.setBase_seq(temp_line.substring(exon.getStart_position() - 1, exon.getEnd_position()));
					}
				}
			}
			file_read.closeFile(genome_fasta);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return out;
	}
	
	void writeFasta(String out_file, ArrayList<String> out) {
		FileWrite file_write = new FileWrite();
		file_write.fileWrite(out_file, out);
	}
	
	/*
	 * exclude_dev means cannot use size, usually equals read length - 1
	 */
	int getExonRandomMax(ArrayList<ExonInfo> exon_list, int exclude_dev, ArrayList<Integer> exon_seq) {
		int out = 0;
		for (int i = 0; i < exon_list.size(); i++) {
			ExonInfo exon = exon_list.get(i);
			int eff_length = exon.getEnd_position() - exon.getStart_position() + 1;
			eff_length -= exclude_dev;
			if (eff_length < 0) {
				eff_length =0;
			}
			out += eff_length;
			exon_seq.add(out);
		}
		return out;
	}
	
	int getExonForwardSpliceRandomMax(ArrayList<ExonInfo> exon_list, int splice_max_length, int splice_min_length, ArrayList<Integer> exon_seq) {
		int out = 0;
		int max_length = splice_max_length;
		int min_length = splice_min_length;
		if (max_length < min_length) {
			max_length = min_length;
			min_length = splice_max_length;
		}
		int forward_num = exon_list.size() - 1 ;
		for (int i = 0; i < forward_num; i++) {
			ExonInfo exon = exon_list.get(i);
			int eff_length = exon.getEnd_position() - exon.getStart_position() + 1;
			if (eff_length >= max_length) {
				eff_length = max_length - min_length + 1;
			}
			else if (eff_length >= min_length) {
				eff_length -= min_length;
				eff_length++;
			}
			else {
				eff_length=0;
			}
			out += eff_length;
			exon_seq.add(out);
		}
		return out;
	}
	
	int[] divideInto(int total, int size) {
		int[] out = new int[size];
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for (int i = 1; i < size; i++) {
			temp.add((int) (Math.random() * total) + 1);
		}
		temp.sort(null);
		out[0] = temp.get(0);
		for (int i = 1; i < size - 1; i++) {
			out[i] = temp.get(i) - temp.get(i - 1);
		}
		out[size - 1] = total - temp.get(size - 2);
		return out;
	}
	
	String getBoundSplice(ArrayList<ExonInfo> exon_list, int exon_num, int offset, int length, StringBuffer note) {
		String out = null;
		ExonInfo exon = exon_list.get(exon_num);
		if (offset < exon.getBase_seq().length()) {
			StringBuffer temp = new StringBuffer();
			temp.append(exon.getBase_seq().substring(offset));
			note.append(exon.getChr_symbol());
			note.append(':');
			note.append(exon.getStart_position() + offset);
			note.append('-');
			note.append(exon.getEnd_position());
			int exon_2 = exon_num;
			int left_length = length - temp.length();
			while (left_length > 0) {
				exon_2 = (int) ((exon_list.size() - exon_2 - 1) * Math.random()) + 1 + exon_2;
				exon = exon_list.get(exon_2);
				String base_seq = exon_list.get(exon_2).getBase_seq();
				note.append('_');
				note.append(exon.getChr_symbol());
				note.append(':');
				note.append(exon.getStart_position());
				note.append('-');
				if (base_seq.length() > left_length) {
					temp.append(base_seq.substring(0, left_length));
					note.append(exon.getStart_position() + left_length - 1);
				}
				else {
					temp.append(base_seq);
					note.append(exon.getEnd_position());
				}
				left_length -= base_seq.length();
			}
			out = temp.toString();
		}
		return out;
	}
	
	String getBoundCirc(ArrayList<ExonInfo> exon_list, int exon_num, int offset, int length, StringBuffer note) {
		String out = null;
		ExonInfo exon = exon_list.get(exon_num);
		if (offset < exon.getBase_seq().length()) {
			StringBuffer temp = new StringBuffer();
			temp.append(exon.getBase_seq().substring(offset));
			note.append("circ:");
			note.append(exon.getChr_symbol());
			note.append(':');
			note.append(exon.getStart_position() + offset);
			note.append('-');
			note.append(exon.getEnd_position());
			int exon_2 = exon_num;
			int left_length = length - temp.length();
			while (left_length > 0) {
				exon_2 = (int) (exon_2 * Math.random());
				exon = exon_list.get(exon_2);
				String base_seq = exon_list.get(exon_2).getBase_seq();
				note.append('_');
				note.append(exon.getChr_symbol());
				note.append(':');
				note.append(exon.getStart_position());
				note.append('-');
				if (base_seq.length() > left_length) {
					temp.append(base_seq.substring(0, left_length));
					note.append(exon.getStart_position() + left_length - 1);
				}
				else {
					temp.append(base_seq);
					note.append(exon.getEnd_position());
				}
				left_length -= base_seq.length();
			}
			out = temp.toString();
		}
		return out;
	}
	
	String getBoundSplicePieces(ArrayList<ExonInfo> exon_list, int exon_num, int offset, int length, StringBuffer note, ArrayList<String> write) {
		String out = null;
		ExonInfo exon = exon_list.get(exon_num);
		if (offset < exon.getBase_seq().length()) {
			StringBuffer temp = new StringBuffer();
			temp.append(exon.getBase_seq().substring(offset));
			note.append(exon.getChr_symbol());
			note.append(':');
			note.append(exon.getStart_position() + offset);
			note.append('-');
			note.append(exon.getEnd_position());
			write.add(note.toString());
			write.add(temp.toString());
			int exon_2 = exon_num;
			int left_length = length - temp.length();
			while (left_length > 0) {
				exon_2 = (int) ((exon_list.size() - exon_2 - 1) * Math.random()) + 1 + exon_2;
				exon = exon_list.get(exon_2);
				String base_seq = exon_list.get(exon_2).getBase_seq();
				note.append('_');
				int note_now = note.length();
				note.append(exon.getChr_symbol());
				note.append(':');
				note.append(exon.getStart_position());
				note.append('-');
				if (base_seq.length() > left_length) {
					temp.append(base_seq.substring(0, left_length));
					note.append(exon.getStart_position() + left_length - 1);
					write.add(">" + note.toString().substring(note_now));
					write.add(base_seq.substring(0, left_length));
				}
				else {
					temp.append(base_seq);
					note.append(exon.getEnd_position());
					write.add(">" + note.toString().substring(note_now));
					write.add(base_seq);
				}
				left_length -= base_seq.length();
			}
			out = temp.toString();
		}
		return out;
	}
	
	int searchMinNoLess(int target, ArrayList<Integer> inc_seq) {
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
	
	static int getTotal_modes() {
		return total_modes;
	}
	
}
