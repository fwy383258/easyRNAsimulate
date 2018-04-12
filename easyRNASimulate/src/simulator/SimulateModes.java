package simulator;

import java.util.ArrayList;
import java.util.HashSet;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class SimulateModes {

	private static int total_modes=4;
	private int scripts_gene=3;
	private int bj_count=20000;
	private int peak_count=2000;
	private double input_dept=0.01;
	
	static void run(ReadInputs args) {
		ArrayList<String> out_reads = null;
		SimulateModes sim = new SimulateModes();
		ArrayList<Chromosome> chr = null;
		int mode = args.getSim_mode();
		switch (mode) {
		case 1:
		case 2:
			ArrayList<ArrayList<Exon>> exons = FileRead.loadFiles(args.getRef_exon_gtf(), args.getRef_genome_file());
			if (mode == 1) {
				out_reads = sim.mode1Uniform(args, exons);
			}
			else {
				out_reads = sim.mode2Binomial(args, exons);
			}
			sim.writeFasta(args.getOut_file(), out_reads);
			break;
		case 3:
			chr = FileRead.loadJunctionFile(args.getRef_bed_file());
			FileRead.loadTranscriptsFile(args, chr);
			out_reads =  sim.mode3TransBin(args, chr);
			sim.writeFasta(args.getOut_file(), out_reads);
			break;
		case 4:
			chr = new ArrayList<>();
			for (int i = 0; i < 25; i++) {
				chr.add(null);
			}
			FileRead.loadTranscriptsFile(args, chr);
			out_reads =  sim.mode4TransExon(args, chr, sim.input_dept);
			break;
		default:
			break;
		}
	}
	
	ArrayList<String> mode1Uniform(ReadInputs args, ArrayList<ArrayList<Exon>> exons){
		ArrayList<String> out = new ArrayList<String>();
		BuildFasta bf = new BuildFasta();
		
		
		int read_count = args.getRead_count();
		int map_rand = (int) (read_count * (1.0 - args.getSplice_input_scale()));
		int circ_rand = read_count - map_rand;
		int splice_rand = (int) (circ_rand * (1.0 - args.getCirc_scale()));
		circ_rand -= splice_rand;
		int[] chr_rand = Method.divideInto(map_rand, exons.size());
		int[] chr_srand = Method.divideInto(splice_rand, exons.size());
		int[] chr_crand = Method.divideInto(splice_rand, exons.size());
		StringBuffer base_seq = new StringBuffer();
		StringBuffer id = new StringBuffer();
		for (int i = 0; i < exons.size(); i++) {
			ArrayList<Exon> exon_list = exons.get(i);
			ArrayList<Integer> exon_seq = new ArrayList<Integer>();
			int random_max = this.getExonRandomMax(exons.get(i), args.getRead_length() - 1, exon_seq);
			System.out.println("Full map " + exon_list.get(0).getChr_symbol() + " : " + exon_seq.size());
			while (chr_rand[i] > 0) {
				int offset = Method.randIntReach(1, random_max);
				int exon_index = this.searchMinNoLess(offset, exon_seq);
				offset += exon_list.get(exon_index).getBase_seq().length() - exon_seq.get(exon_index) - args.getRead_length();
				base_seq.setLength(0);
				id.setLength(0);
				id.append('>');
				int	left_length = bf.buildFirstSplice(exon_list.get(exon_index), offset, args.getRead_length(), base_seq, id);
				if (left_length > 0) {
					System.out.println("Error: Full map left");
				}
				out.add(id.toString());
				out.add(base_seq.toString());
				chr_rand[i] --;
			}
			exon_seq = new ArrayList<Integer>();
			random_max = this.getExonSpliceRandomMax(exons.get(i), args.getMax_splice(), args.getMin_splice(), exon_seq);
			System.out.println("Splice " + exon_list.get(0).getChr_symbol() + " : " + exon_seq.size());
			while (chr_srand[i] > 0) {
				int offset = Method.randIntReach(1, random_max);
				int exon_1 = this.searchMinNoLess(offset, exon_seq);
				offset += exon_list.get(exon_1).getBase_seq().length() - exon_seq.get(exon_1) - args.getMin_splice();
				base_seq.setLength(0);
				id.setLength(0);
				id.append('>');
				int left_length = bf.buildFirstSplice(exon_list.get(exon_1), offset, args.getRead_length(), base_seq, id);
				bf.addAllBoundForwardSplice(base_seq, id, left_length, exon_1, exon_list);
				out.add(id.toString());
				out.add(base_seq.toString());
				chr_srand[i] --;
			}
			while (chr_crand[i] > 0) {
				int offset = Method.randIntReach(1, random_max);
				int exon_1 = this.searchMinNoLess(offset, exon_seq);
				offset += exon_list.get(exon_1).getBase_seq().length() - exon_seq.get(exon_1) - args.getMin_splice();
				base_seq.setLength(0);
				id.setLength(0);
				id.append(">circ:");
				int left_length = bf.buildFirstSplice(exon_list.get(exon_1), offset, args.getRead_length(), base_seq, id);
				int exon_2 = Method.randIntReach(0, exon_1);
				left_length = bf.addBoundSplice(base_seq, id, left_length, exon_list.get(exon_2));
				bf.addAllBoundForwardSplice(base_seq, id, left_length, exon_2, exon_list);
				out.add(id.toString());
				out.add(base_seq.toString());
				chr_crand[i] --;
			}
		}
		return out;
	}
	
	ArrayList<String> mode2Binomial(ReadInputs args, ArrayList<ArrayList<Exon>> exons){
		ArrayList<String> out = null;
		BuildFasta bf = new BuildFasta();
		
		int peak_reads = (int) (args.getRead_count() * (Math.random() * 0.2 + 0.1)); 
		args.setRead_count(args.getRead_count() - peak_reads);
		out = this.mode1Uniform(args, exons);
		int peak_num = args.getPeak_num();
		int peak_ava = peak_reads / peak_num / 2;
		peak_reads += - peak_num * peak_ava;
		int[] peak_count = Method.divideInto(peak_reads, peak_num);
		int[] peak_rand = Method.divideInto(peak_num, exons.size());
		
		StringBuffer base_seq = new StringBuffer();
		StringBuffer id = new StringBuffer();
		peak_num=0;
		for (int chr_num = 0; chr_num < exons.size(); chr_num ++) {
			ArrayList<Exon> exon_list = exons.get(chr_num);
			ArrayList<Integer> exon_seq = new ArrayList<>();
			int full_max = this.getExonRandomMax(exon_list, args.getRead_length() - 1, exon_seq);
			while (peak_rand[chr_num] > 0) {
				peak_count[peak_num] += peak_ava;
				int rand_mode = Method.randIntReach(1, 3);
				if (rand_mode == 1) {
					int exon_index = this.searchMinNoLess(Method.randIntReach(1, full_max), exon_seq);
					int length = exon_list.get(exon_index).getBase_seq().length() - args.getRead_length();
					int left_point = 0;
					BinomialDistribution bd = null;
					if (length > args.getRead_length()) {
						left_point = Method.randIntReach(0, length - args.getRead_length() + 1);
						bd = new BinomialDistribution(args.getRead_length(), 0.5);
					}
					else {
						bd = new BinomialDistribution(length, 0.5);
					}
					while(peak_count[peak_num] > 0) {
						int offset = left_point + bd.sample();
						base_seq.setLength(0);
						id.setLength(0);
						id.append(">peak");
						id.append(peak_count[peak_num]);
						id.append(':');
						int	left_length = bf.buildFirstSplice(exon_list.get(exon_index), offset, args.getRead_length(), base_seq, id);
						if (left_length > 0) {
							System.out.println("Error: Full map left peak");
						}
						out.add(id.toString());
						out.add(base_seq.toString());
						peak_count[peak_num]--;
					}
				}
				else if(rand_mode == 2) {
					int exon_1 = this.searchMinNoLess(Method.randIntReach(1, full_max), exon_seq);
					int left_point = exon_list.get(exon_1).getBase_seq().length() - args.getRead_length() + 1;
					BinomialDistribution bd = new BinomialDistribution(args.getRead_length(), 0.5);
					int exon_2 = exon_seq.size() - 1;
					if (exon_1 < exon_2) {
						exon_2 = Method.randIntReach(exon_1 + 1, exon_2);
						while(peak_count[peak_num] > 0) {
							int offset = left_point + bd.sample();
							base_seq.setLength(0);
							id.setLength(0);
							id.append(">peak");
							id.append(peak_count[peak_num]);
							id.append(':');
							int	left_length = bf.buildFirstSplice(exon_list.get(exon_1), offset, args.getRead_length(), base_seq, id);
							left_length = bf.addBoundSplice(base_seq, id, left_length, exon_list.get(exon_2));
							bf.addAllBoundForwardSplice(base_seq, id, left_length, exon_2, exon_list);
							out.add(id.toString());
							out.add(base_seq.toString());
							peak_count[peak_num]--;
						}
					}
					else {
						while(peak_count[peak_num] > 0) {
							int offset = left_point + bd.sample();
							base_seq.setLength(0);
							id.setLength(0);
							id.append(">peak");
							id.append(peak_count[peak_num]);
							id.append(":circ:");
							int	left_length = bf.buildFirstSplice(exon_list.get(exon_1), offset, args.getRead_length(), base_seq, id);
							left_length = bf.addBoundSplice(base_seq, id, left_length, exon_list.get(exon_2));
							bf.addAllBoundForwardSplice(base_seq, id, left_length, exon_2, exon_list);
							out.add(id.toString());
							out.add(base_seq.toString());
							peak_count[peak_num]--;
						}
					}
				}
				else {
					int exon_1 = this.searchMinNoLess(Method.randIntReach(1, full_max), exon_seq);
					int left_point = exon_list.get(exon_1).getBase_seq().length() - args.getRead_length() + 1;
					BinomialDistribution bd = new BinomialDistribution(args.getRead_length(), 0.5);
					int exon_2 =Method.randIntReach(0, exon_1);
					while(peak_count[peak_num] > 0) {
						int offset = left_point + bd.sample();
						base_seq.setLength(0);
						id.setLength(0);
						id.append(">peak");
						id.append(peak_count[peak_num]);
						id.append(":circ:");
						int	left_length = bf.buildFirstSplice(exon_list.get(exon_1), offset, args.getRead_length(), base_seq, id);
						left_length = bf.addBoundSplice(base_seq, id, left_length, exon_list.get(exon_2));
						bf.addAllBoundForwardSplice(base_seq, id, left_length, exon_2, exon_list);
						out.add(id.toString());
						out.add(base_seq.toString());
						peak_count[peak_num]--;
					}
				}
				peak_rand[chr_num]--;
				peak_num++;
			}
		}
		
		return out;
	}
	
	ArrayList<String> mode3TransBin(ReadInputs args, ArrayList<Chromosome> chr_list){
		ArrayList<String> out = new ArrayList<>(64);
		ArrayList<String> out_ip = new ArrayList<>(64);
		HashSet<Transcript> temp = new HashSet<>();
		BuildFasta bf = new BuildFasta();
		int warn_count = 0;
		
		for(int i=0; i < chr_list.size(); i++) {
			Chromosome chr = chr_list.get(i);
			if (chr != null) {
				for (int j=0; j < chr.getGenes().size(); j++) {
					Gene gene = chr.getGene(j);
					if (gene.getScripts().size() < this.scripts_gene) {
						temp.addAll(gene.getScripts());
					}
					else {
						Method.getNoneRepeat(gene.getScripts(), this.scripts_gene, temp);
					}
					for (Transcript the_script : temp) {
						PoissonDistribution p = new PoissonDistribution(5.30);
						bf.getScriptReads(out, out_ip, null, the_script, p.sample(), 0, args.getRead_length(), false, false);
					}
					temp.clear();
				}

				int bj_count = chr.getAll_junctions().size() / 100;
				while (bj_count > 0) {
					int bj_index = Method.getRandElement(chr.getAll_junctions());
					BackJunction bj = Method.getRandElement(chr.getJunction_map().get(bj_index));
					if (bj.getBoth_scripts().size() == 0) {
						warn_count++;
						continue;
					}
					warn_count--;
					Transcript script = Method.getRandElement(bj.getBoth_scripts());
					Transcript the_script = new Transcript(null, script.getId(), new ArrayList<>(), 0, 0, 2 * args.getRead_length() + 2, false);
					Exon exon = new Exon();
					exon.setBase_seq(bj.getEnd_seq());
					exon.setChr_symbol(chr.getId());
					exon.setStart(bj.getEnd() - args.getRead_length());
					exon.setEnd(bj.getEnd());
					the_script.addExon(exon);
					exon = new Exon();
					exon.setBase_seq(bj.getStart_seq());
					exon.setChr_symbol(chr.getId());
					exon.setStart(bj.getStart());
					exon.setEnd(bj.getStart() + args.getRead_length());
					the_script.addExon(exon);
					the_script.setCirc_flag(true);
					PoissonDistribution p = new PoissonDistribution(3.50);
					bf.getScriptReads(out, out_ip, null, the_script, p.sample(), 0, args.getRead_length(), false, false);
					bj_count--;
				}
			}
		}
		System.out.println("Total warn in back junction: " + warn_count);
		return out;
	}
	
	ArrayList<String> mode4TransExon(ReadInputs args, ArrayList<Chromosome> chr_list, double read_depth){
		ArrayList<String> out = new ArrayList<>();
		ArrayList<String> out_ip = new ArrayList<>();
		HashSet<Transcript> peak_scripts = new HashSet<>();
		HashSet<Transcript> bj_scripts = this.getRandBJ(chr_list, this.bj_count, args.getRead_length() * 2, peak_scripts, this.peak_count);
		HashSet<Transcript> temp = new HashSet<>();
		BuildFasta bf = new BuildFasta();
		FileWrite fw = new FileWrite();
		String prefix = args.getOut_file().substring(0, args.getOut_file().lastIndexOf('.'));
		String ip_file = prefix + "_IP.fa";
		int suc_peak = 0;
		ArrayList<String> peak_out = new ArrayList<>();
		peak_out.add("Chr\tStart\tEnd\tID");
		fw.fileWrite(args.getOut_file(), out);
		fw.fileWrite(ip_file, out);
		
		for(int i=0; i < chr_list.size(); i++) {
			Chromosome chr = chr_list.get(i);
			if (chr != null) {
				for (int j=0; j < chr.getGenes().size(); j++) {
					Gene gene = chr.getGene(j);
					ArrayList<Transcript> temp_list = new ArrayList<>();
					temp_list.addAll(gene.getScripts());
					temp_list.removeAll(bj_scripts);
					temp_list.removeAll(peak_scripts);
					if (temp_list.size() > this.scripts_gene) {
						Method.getNoneRepeat(temp_list, this.scripts_gene, temp);
					}
					else {
						temp.addAll(temp_list);
					}
					for (Transcript the_script : temp) {
						PoissonDistribution p = new PoissonDistribution(read_depth * the_script.getBase_sum());
						PoissonDistribution p_ip = new PoissonDistribution(0.5 * read_depth * the_script.getBase_sum());
						bf.getScriptReads(out, out_ip, null,the_script, p.sample(), p_ip.sample(), args.getRead_length(), false, false);
					}
					temp.clear();
				}
			}
			fw.fileAppend(args.getOut_file(), out);
			fw.fileAppend(ip_file, out_ip);
			out.clear();
			out_ip.clear();
		}
		int suc_bj = 0;
		ArrayList<String> circ_out = new ArrayList<>();
		for (Transcript the_script : bj_scripts) {
			PoissonDistribution p = new PoissonDistribution(1.3 * read_depth * the_script.getBase_sum());;
			PoissonDistribution p_ip = null; 
			int p_count = 0;
			int ip_count = 0;
			boolean peak_flag = peak_scripts.contains(the_script);
			if (peak_flag) {
				p_ip = new PoissonDistribution(130.0 * read_depth * the_script.getBase_sum());
				while (p_count >= ip_count) {
					p_count = p.sample();
					ip_count = p_ip.sample();
				}
			}
			else {
				p_ip = new PoissonDistribution(0.65 * read_depth * the_script.getBase_sum());
				p_count = p.sample();
				ip_count = p_ip.sample();
			}
			if(bf.getScriptReads(out, out_ip, peak_out, the_script, p_count, ip_count, args.getRead_length(), true, peak_flag)) {
				if (peak_flag) {
					suc_peak++;
					peak_scripts.remove(the_script);
				}
				p_count = p.sample();
				ip_count = p_ip.sample();
				bf.getScriptReads(out, out_ip, null, the_script, p_count, ip_count, args.getRead_length(), false, false);
				circ_out.add(the_script.getId() + "\t" + the_script.getExon(0).getChr_symbol() + "\t" + the_script.getExon(1).getStart() + "\t" + the_script.getExon(0).getEnd());
				suc_bj++;
			}
			fw.fileAppend(args.getOut_file(), out);
			fw.fileAppend(ip_file, out_ip);
			out.clear();
			out_ip.clear();
		}
		for (Transcript the_script : peak_scripts) {
			PoissonDistribution p = new PoissonDistribution(read_depth * the_script.getBase_sum());
			PoissonDistribution p_ip = new PoissonDistribution(100.0 * read_depth * the_script.getBase_sum());
			int p_count = 0;
			int ip_count = 0;
			while (p_count >= ip_count) {
				p_count = p.sample();
				ip_count = p_ip.sample();
			}
			if(bf.getScriptReads(out, out_ip, peak_out, the_script, p_count, ip_count, args.getRead_length(), false, true)) {
				suc_peak++;
				fw.fileAppend(args.getOut_file(), out);
				fw.fileAppend(ip_file, out_ip);
				out.clear();
				out_ip.clear();
			}
		}
		System.out.println("Circ succeed : " + suc_bj);
		System.out.println("Peak succeed : " + suc_peak);
		fw.fileWrite(prefix + ".circ", circ_out);
		fw.fileWrite(prefix + ".peak", peak_out);
		return out;
	}
	
	void writeFasta(String out_file, ArrayList<String> out) {
		FileWrite file_write = new FileWrite();
		file_write.fileWrite(out_file, out);
	}
	
	/*
	 * exclude_dev means cannot use size, usually equals read length - 1
	 */
	int getExonRandomMax(ArrayList<Exon> exon_list, int exclude_dev, ArrayList<Integer> exon_seq) {
		int out = 0;
		for (int i = 0; i < exon_list.size(); i++) {
			Exon exon = exon_list.get(i);
			int eff_length = exon.getEnd() - exon.getStart() + 1;
			eff_length -= exclude_dev;
			if (eff_length < 0) {
				eff_length =0;
			}
			out += eff_length;
			exon_seq.add(out);
		}
		return out;
	}
	
	int getExonSpliceRandomMax(ArrayList<Exon> exon_list, int splice_max_length, int splice_min_length, ArrayList<Integer> exon_seq) {
		int out = 0;
		int max_length = splice_max_length;
		int min_length = splice_min_length;
		if (max_length < min_length) {
			max_length = min_length;
			min_length = splice_max_length;
		}
		int forward_num = exon_list.size() - 1 ;
		for (int i = 0; i < forward_num; i++) {
			Exon exon = exon_list.get(i);
			int eff_length = exon.getEnd() - exon.getStart() + 1;
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
	
	String getBoundSplice(ArrayList<Exon> exon_list, int exon_num, int offset, int length, StringBuffer note) {
		String out = null;
		Exon exon = exon_list.get(exon_num);
		if (offset < exon.getBase_seq().length()) {
			StringBuffer temp = new StringBuffer();
			temp.append(exon.getBase_seq().substring(offset));
			note.append(exon.getChr_symbol());
			note.append(':');
			note.append(exon.getStart() + offset);
			note.append('-');
			note.append(exon.getEnd());
			int exon_2 = exon_num;
			int left_length = length - temp.length();
			while (left_length > 0) {
				exon_2 = (int) ((exon_list.size() - exon_2 - 1) * Math.random()) + 1 + exon_2;
				exon = exon_list.get(exon_2);
				String base_seq = exon_list.get(exon_2).getBase_seq();
				note.append('_');
				note.append(exon.getChr_symbol());
				note.append(':');
				note.append(exon.getStart());
				note.append('-');
				if (base_seq.length() > left_length) {
					temp.append(base_seq.substring(0, left_length));
					note.append(exon.getStart() + left_length - 1);
				}
				else {
					temp.append(base_seq);
					note.append(exon.getEnd());
				}
				left_length -= base_seq.length();
			}
			out = temp.toString();
		}
		return out;
	}
	
	String getBoundCirc(ArrayList<Exon> exon_list, int exon_num, int offset, int length, StringBuffer note) {
		String out = null;
		Exon exon = exon_list.get(exon_num);
		if (offset < exon.getBase_seq().length()) {
			StringBuffer temp = new StringBuffer();
			temp.append(exon.getBase_seq().substring(offset));
			note.append("circ:");
			note.append(exon.getChr_symbol());
			note.append(':');
			note.append(exon.getStart() + offset);
			note.append('-');
			note.append(exon.getEnd());
			int exon_2 = exon_num;
			int left_length = length - temp.length();
			while (left_length > 0) {
				exon_2 = (int) (exon_2 * Math.random());
				exon = exon_list.get(exon_2);
				String base_seq = exon_list.get(exon_2).getBase_seq();
				note.append('_');
				note.append(exon.getChr_symbol());
				note.append(':');
				note.append(exon.getStart());
				note.append('-');
				if (base_seq.length() > left_length) {
					temp.append(base_seq.substring(0, left_length));
					note.append(exon.getStart() + left_length - 1);
				}
				else {
					temp.append(base_seq);
					note.append(exon.getEnd());
				}
				left_length -= base_seq.length();
			}
			out = temp.toString();
		}
		return out;
	}
	
	String getBoundSplicePieces(ArrayList<Exon> exon_list, int exon_num, int offset, int length, StringBuffer note, ArrayList<String> write) {
		String out = null;
		Exon exon = exon_list.get(exon_num);
		if (offset < exon.getBase_seq().length()) {
			StringBuffer temp = new StringBuffer();
			temp.append(exon.getBase_seq().substring(offset));
			note.append(exon.getChr_symbol());
			note.append(':');
			note.append(exon.getStart() + offset);
			note.append('-');
			note.append(exon.getEnd());
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
				note.append(exon.getStart());
				note.append('-');
				if (base_seq.length() > left_length) {
					temp.append(base_seq.substring(0, left_length));
					note.append(exon.getStart() + left_length - 1);
					write.add(">" + note.toString().substring(note_now));
					write.add(base_seq.substring(0, left_length));
				}
				else {
					temp.append(base_seq);
					note.append(exon.getEnd());
					write.add(">" + note.toString().substring(note_now));
					write.add(base_seq);
				}
				left_length -= base_seq.length();
			}
			out = temp.toString();
		}
		return out;
	}
	
	HashSet<Transcript> getRandBJ(ArrayList<Chromosome> chr_list, int count, int limit, HashSet<Transcript> peak_scripts, int peak_count){
		HashSet<Transcript> out = new HashSet<>();
		ArrayList<Transcript> temp = new ArrayList<>();
		for (int i=0; i < chr_list.size(); i++) {
			if (chr_list.get(i) != null) {
				for (int j=0; j < chr_list.get(i).getGenes().size(); j++) {
					Gene the_gene = chr_list.get(i).getGene(j);
					temp.addAll(the_gene.getScripts());
				}
			}
		}
		if (temp.size() > count) {
			Method.getNoneRepeat(temp, count, out);
			temp.removeAll(out);
			ArrayList<Transcript> circ = new ArrayList<>();
			circ.addAll(out);
			int peak_circ = peak_count * count / temp.size();
			Method.getNoneRepeat(circ, peak_circ, peak_scripts);
			while(temp.size() > 0 && peak_scripts.size() < peak_count) {
				Transcript the_script = temp.get(Method.randIntReach(0, temp.size() - 1));
				if (the_script.getBase_sum() > limit) {
					peak_scripts.add(the_script);
				}
				temp.remove(the_script);
			}
			if (peak_scripts.size() < peak_count) {
				System.out.println("Warning: Not enough elements for picking");
			}
		}
		else {
			out.addAll(temp);
			System.out.println("Warning: wrong BJ count, all scripts have BJ");
		}
		for (Transcript script : out) {
			int bj_index = Method.randIntReach(0, script.getExons().size() - 1);
			Exon bj_exon = script.getExon(bj_index);
			int length = 0;
			ArrayList<Exon> replace_list = new ArrayList<>();
			replace_list.add(bj_exon);
			length += bj_exon.getBase_seq().length();
			for (int i=0; i < bj_index; i++) {
				if (Math.random() < 0.667 && replace_list.size() < 4) {
					length += script.getExon(i).getBase_seq().length();
					replace_list.add(script.getExon(i));
				}
			}
			while (length < limit || replace_list.size() < 2) {
				length += bj_exon.getBase_seq().length();
				replace_list.add(bj_exon);
			}
			script.setBase_sum(length);
			script.setExons(replace_list);
			script.setCirc_flag(true);
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
