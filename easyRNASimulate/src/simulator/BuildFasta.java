package simulator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class BuildFasta {
	
	private int fragment_length = 200;
	
	boolean getScriptReads(ArrayList<String> out, ArrayList<String> out_ip, ArrayList<String> circ_report, ArrayList<String> peak_report, Transcript script, int count, int count_ip, int length, boolean circ_flag, boolean peak_flag) {
		StringBuffer id = new StringBuffer();
		StringBuffer base_seq = new StringBuffer();
		
		int base_sum = script.getBase_sum();
		if (base_sum < length) {
			return false;
		}
		id.append(">");
		int start_index = 0;
		int end_index = base_sum - length;
		if (script.isCirc_flag()) {
			int bj_index = script.getExon(0).getEnd() - script.getExon(0).getStart();
			if (circ_flag) {
				start_index = bj_index - length + 1;
				end_index = bj_index;
				if (start_index < 0) {
					end_index -= start_index;
					start_index = 0;
				}
				if (end_index > base_sum - length) {
					start_index -= end_index - base_sum + length;
					end_index = base_sum - length;
				}
				if (start_index < 0) {
					return false;
				}
				if (peak_flag) {
					id.append("peak:");
				}
				id.append("circ:");
			}
			else {
				int front_l = bj_index - 2 * length + 2;
				if (front_l < 0) {
					front_l = 0;
				}
				int back_l = end_index - length - bj_index;
				if (back_l < 0) {
					back_l = 0;
				}
				if (front_l + back_l == 0) {
					return false;
				}
				int rand = Method.randIntReach(0, front_l + back_l - 1);
				if (rand >= front_l) {
					rand += bj_index + 1 - front_l;
				}
				start_index = rand;
				end_index = rand + length - 1;
				if (peak_flag) {
					if (start_index - length / 4 >= 0 && end_index + length / 4 < base_sum - length) {
						rand = Method.randIntReach(- length/4, length/4);
						start_index += rand;
						end_index += rand;
					}
					id.append("peak:");
				}
				id.append("circn:");
			}
		}
		else if (peak_flag) {
			if (end_index - length > start_index) {
				int rand = Method.randIntReach(start_index, end_index - length);
				start_index = rand;
				end_index = rand + length;
			}
			id.append("peak:");
		}
		id.append(script.getId());
		id.append(':');
		int id_fix = id.length();
		Bed12 report_line = null;
		int len = 0;
		if (circ_report != null && circ_flag && script.isCirc_flag()) {
			report_line = new Bed12();
			report_line.setName(script.getId());
			report_line.setChr(script.getExon(0).getChr_symbol());
			report_line.setStrand(script.getStrand());
			report_line.setStart(script.getExon(1).getStart());
			report_line.setEnd(script.getExon(0).getEnd());
			for (int i = 0; i < script.getExons().size(); i++) {
				Exon exon = script.getExon(i);
				if (start_index < exon.getEnd() - exon.getStart() + 1 + len) {
					report_line.setBlock_count(report_line.getBlock_count() + 1);
					if (start_index >= len) {
						report_line.getBlock_starts().add(exon.getStart() + start_index - report_line.getStart());
						report_line.getBlock_sizes().add(exon.getEnd() - exon.getStart() + 1 + len - start_index);
					}
					else {
						report_line.getBlock_starts().add(exon.getStart() - report_line.getStart());
						report_line.getBlock_sizes().add(exon.getEnd() - exon.getStart() + 1);
					}
					if (end_index + length <= exon.getEnd() - exon.getStart() + 1 + len) {
						report_line.getBlock_sizes().set(report_line.getBlock_count() - 1, end_index + length - len);
						report_line.setThick_start(report_line.getStart());
						report_line.setThick_end(report_line.getEnd());
						break;
					}
				}
				len += exon.getEnd() - exon.getStart() + 1;
			}
			circ_report.add(report_line.toString());
		}
		if (peak_report != null && peak_flag) {
			if (report_line == null) {
				report_line = new Bed12();
				report_line.setName(script.getId());
				report_line.setChr(script.getExon(0).getChr_symbol());
				report_line.setStrand(script.getStrand());
				if (script.isCirc_flag() && circ_flag) {
					report_line.setStart(script.getExon(1).getStart());
					report_line.setEnd(script.getExon(0).getEnd());
				}
				for (int i = 0; i < script.getExons().size(); i++) {
					Exon exon = script.getExon(i);
					if (start_index < exon.getEnd() - exon.getStart() + 1 + len) {
						report_line.setBlock_count(report_line.getBlock_count() + 1);
						if (start_index >= len) {
							if (report_line.getStart() == 0) {
								report_line.setStart(exon.getStart() + start_index - len);
								report_line.getBlock_starts().add(0);
								report_line.getBlock_sizes().add(exon.getEnd() - report_line.getStart() + 1);
							}
							else {
								report_line.getBlock_starts().add(exon.getStart() + start_index - report_line.getStart());
								report_line.getBlock_sizes().add(exon.getEnd() - exon.getStart() + 1 + len - start_index);
							}
						}
						else {
							report_line.getBlock_starts().add(exon.getStart() - report_line.getStart());
							report_line.getBlock_sizes().add(exon.getEnd() - exon.getStart() + 1);
						}
						if (end_index + length <= exon.getEnd() - exon.getStart() + 1 + len) {
							if (report_line.getEnd() == 0) {
								report_line.setEnd(exon.getStart() + end_index + length - len);
							}
							report_line.getBlock_sizes().set(report_line.getBlock_count() - 1, end_index + length - len);
							report_line.setThick_start(report_line.getStart());
							report_line.setThick_end(report_line.getEnd());
							break;
						}
					}
					len += exon.getEnd() - exon.getStart() + 1;
				}
			}
			peak_report.add(report_line.toString());
		}
		this.getScriptReads(out, script, count, id, id_fix, base_seq, start_index, end_index, length);
		this.getScriptReads(out_ip, script, count_ip, id, id_fix, base_seq, start_index, end_index, length);
//		HashMap<Integer, Integer> count_id = new HashMap<Integer, Integer>();
//		while (count > 0) {
//			id.setLength(id_fix);
//			base_seq.setLength(0);
//			int rand = Method.randIntReach(start_index, end_index);
//			if (count_id.containsKey(rand)) {
//				count_id.put(rand, count_id.get(rand) + 1);
//			}
//			else {
//				count_id.put(rand, 1);
//			}
//			id.append(count_id.get(rand));
//			id.append(':');
//			for (int i=0; i < script.getExons().size(); i++) {
//				Exon exon = script.getExon(i);
//				int l = exon.getEnd() - exon.getStart() + 1;
//				if (rand < l) {
//					int left_length = this.buildFirstSplice(exon, rand, length, base_seq, id);
//					int j = i;
//					while (left_length > 0) {
//						j++;
//						left_length = this.addBoundSplice(base_seq, id, left_length, script.getExon(j));
//					}
//					break;
//				}
//				else {
//					rand -= l;
//				}
//			}			
//			out.add(id.toString());
//			out.add(base_seq.toString());
//			count--;
//		}
//		while (count_ip > 0) {
//			id.setLength(id_fix);
//			base_seq.setLength(0);
//			int rand = Method.randIntReach(start_index, end_index);
//			if (count_id.containsKey(rand)) {
//				count_id.put(rand, count_id.get(rand) + 1);
//			}
//			else {
//				count_id.put(rand, 1);
//			}
//			id.append(count_id.get(rand));
//			id.append(':');
//			for (int i=0; i < script.getExons().size(); i++) {
//				Exon exon = script.getExon(i);
//				int l = exon.getEnd() - exon.getStart() + 1;
//				if (rand < l) {
//					int left_length = this.buildFirstSplice(exon, rand, length, base_seq, id);
//					int j = i;
//					while (left_length > 0) {
//						j++;
//						left_length = this.addBoundSplice(base_seq, id, left_length, script.getExon(j));
//					}
//					break;
//				}
//				else {
//					rand -= l;
//				}
//			}			
//			out_ip.add(id.toString());
//			out_ip.add(base_seq.toString());
//			count_ip--;
//		}
		return true;
	}
	
	public void getForwardReads(ArrayList<Read> out, ArrayList<Read> out_ip, Transcript script, int count, int count_ip, int length) {
		int start_index = 0;
		int end_index = script.getBase_sum() - length;
		for (int i = 0; i < count; i++) {
			int pos = Method.randIntReach(start_index, end_index);
			Read the_read = new Read();
			for (int j = 0; j < script.getExons().size(); j++) {
				Exon exon = script.getExon(j);
				int l = exon.getEnd() - exon.getStart() + 1;
				if (pos < l) {
					int left_length = length - l + pos;
					the_read.getPositions().add(exon.getStart() + pos);
					while (left_length > 0) {
						the_read.getPositions().add(exon.getEnd());
						j++;
						exon = script.getExon(j);
						the_read.getPositions().add(exon.getStart());
						left_length -= exon.getEnd() - exon.getStart() + 1;
					}
					the_read.getPositions().add(left_length + exon.getEnd());
					break;
				}
				else {
					pos -= l;
				}
			}
			the_read.setId(script.getId());
			out.add(the_read);
		}
		
		for (int i = 0; i < count_ip; i++) {
			int pos = Method.randIntReach(start_index, end_index);
			Read the_read = new Read();
			for (int j = 0; j < script.getExons().size(); j++) {
				Exon exon = script.getExon(j);
				int l = exon.getEnd() - exon.getStart() + 1;
				if (pos < l) {
					int left_length = length - l + pos;
					the_read.getPositions().add(exon.getStart() + pos);
					while (left_length > 0) {
						the_read.getPositions().add(exon.getEnd());
						j++;
						exon = script.getExon(j);
						the_read.getPositions().add(exon.getStart());
						left_length -= exon.getEnd() - exon.getStart() + 1;
					}
					the_read.getPositions().add(left_length + exon.getEnd());
					break;
				}
				else {
					pos -= l;
				}
			}
			the_read.setId(script.getId());
			out_ip.add(the_read);
		}
	}
	
	public Bed12 getCircReads(ArrayList<Read> out, ArrayList<Read> out_ip, Transcript script, int count, int count_ip, int length) {
		if (script.isCirc_flag()) {
			boolean intro_flag = false;
			int start_index = 0;
			int end_index = 0;
			start_index = Method.randIntReach(0, script.getExons().size() - 1);
			end_index = Method.randIntReach(0, script.getExons().size() - 1);
			if (start_index == end_index) {
				Exon exon = script.getExon(start_index);
				if (exon.getEnd() - exon.getStart() >= 199) {
					script.getExons().clear();
					script.addExon(exon);
					script.setBase_sum(exon.getEnd() - exon.getStart() + 1);
				}
				else {
					int rand = Method.randIntReach(0, script.getExons().size() - 2);
					if (rand >= start_index) {
						end_index = rand + 1;
						
					}
					else {
						start_index = rand;
					}
				}
			}
			else if (start_index > end_index){
				start_index ^= end_index;
				end_index ^= start_index;
				start_index ^= end_index;
			}
			
			if (start_index < end_index) {
				ArrayList<Exon> exons = new ArrayList<>();
				int base_sum = 0;
				int last_index = start_index;
				exons.add(script.getExon(end_index));
				base_sum += script.getExon(end_index).getEnd() - script.getExon(end_index).getStart() + 1;
				exons.add(script.getExon(start_index));
				base_sum += script.getExon(start_index).getEnd() - script.getExon(start_index).getStart() + 1;
				for (int i = start_index + 1; i < end_index; i++) {
					if (Math.random() < 0.5) {
						exons.add(script.getExon(i));
						base_sum += script.getExon(i).getEnd() - script.getExon(i).getStart() + 1;
						last_index = i;
					}
				}
				if (last_index == end_index - 1) {
					if (base_sum < 200 || Math.random() < 0.5) {
						base_sum += script.getExon(end_index).getStart() - script.getExon(last_index).getEnd() - 1;
						script.getExon(last_index).setEnd(script.getExon(end_index).getStart() - 1);
						intro_flag = true;
					}
				}
				script.setExons(exons);
				script.setBase_sum(base_sum);
			}
			
			start_index = 0;
			end_index = script.getBase_sum() - 1;
			for (int i = 0; i < count; i++) {
				int pos = Method.randIntReach(start_index, end_index);
				Read the_read = new Read();
				for (int j = 0; ; j++) {
					Exon exon = script.getExon(j);
					int l = exon.getEnd() - exon.getStart() + 1;
					if (pos < l) {
						int left_length = length - l + pos;
						the_read.getPositions().add(exon.getStart() + pos);
						while (left_length > 0) {
							j++;
							j %= script.getExons().size();
							if (exon.getEnd() + 1 == script.getExon(j).getStart()) {
								exon = script.getExon(j);
							}
							else {
								the_read.getPositions().add(exon.getEnd());
								exon = script.getExon(j);
								the_read.getPositions().add(exon.getStart());
							}
							left_length -= exon.getEnd() - exon.getStart() + 1;
						}
						the_read.getPositions().add(left_length + exon.getEnd());
						break;
					}
					else {
						pos -= l;
					}
				}
				the_read.setId(script.getId());
				the_read.setCirc_flag(true);
				out.add(the_read);
			}
			
			for (int i = 0; i < count_ip; i++) {
				int pos = Method.randIntReach(start_index, end_index);
				Read the_read = new Read();
				for (int j = 0; ; j++) {
					Exon exon = script.getExon(j);
					int l = exon.getEnd() - exon.getStart() + 1;
					if (pos < l) {
						int left_length = length - l + pos;
						the_read.getPositions().add(exon.getStart() + pos);
						while (left_length > 0) {
							j++;
							j %= script.getExons().size();
							if (exon.getEnd() + 1 == script.getExon(j).getStart()) {
								exon = script.getExon(j);
							}
							else {
								the_read.getPositions().add(exon.getEnd());
								exon = script.getExon(j);
								the_read.getPositions().add(exon.getStart());
							}
							left_length -= exon.getEnd() - exon.getStart() + 1;
						}
						the_read.getPositions().add(left_length + exon.getEnd());
						break;
					}
					else {
						pos -= l;
					}
				}
				the_read.setId(script.getId());
				the_read.setCirc_flag(true);
				out_ip.add(the_read);
			}
			
			
			Bed12 record = new Bed12();
			record.setName(script.getId());
			record.setChr(script.getExon(0).getChr_symbol());
			record.setStrand(script.getStrand());
			if (script.getExons().size() > 1) {
				record.setStart(script.getExon(1).getStart());
			}
			else {
				record.setStart(script.getExon(0).getStart());
			}
			record.setEnd(script.getExon(0).getEnd());
			record.setBlock_count(script.getExons().size());
			record.setThick_start(record.getStart());
			record.setThick_end(record.getEnd());
			for (int i = 1; i < script.getExons().size(); i++) {
				Exon exon = script.getExon(i);
				record.getBlock_sizes().add(exon.getEnd() - exon.getStart() + 1);
				record.getBlock_starts().add(exon.getStart() - record.getStart());
			}
			if (!intro_flag) {
				record.getBlock_sizes().add(script.getExon(0).getEnd() - script.getExon(0).getStart() + 1);
				record.getBlock_starts().add(script.getExon(0).getStart() - record.getStart());
			}
			else {
				record.getBlock_sizes().set(record.getBlock_sizes().size() - 1, script.getExon(0).getEnd() - script.getExon(script.getExons().size() -1).getStart() + 1);
				record.setBlock_count(record.getBlock_count() - 1);
			}
			return record;
		}
		else {
			System.out.println("Failed to add circ in " + script.getId());
			return null;
		}
	}
	
	public Bed12 getPeakReads(ArrayList<Read> out, ArrayList<Read> out_ip, Transcript script, int count,  int count_ip, int length, boolean circ_peak) {
		int start_index = 0;
		int end_index = 0;
		if (script.isCirc_flag()) {
			end_index = script.getBase_sum() - 1;
			if (circ_peak) {
				start_index = script.getExon(0).getEnd() - script.getExon(0).getStart() - fragment_length;
			}
			else {
				start_index = Method.randIntReach(start_index, end_index);
			}
			end_index = start_index + fragment_length;
			for (int i = 0; i < count; i++) {
				int pos = Method.randIntReach(start_index, end_index);
				if (pos < 0) {
					pos += script.getBase_sum();
				}
				Read the_read = new Read();
				for (int j = 0; ; j++) {
					j %= script.getExons().size();
					Exon exon = script.getExon(j);
					int l = exon.getEnd() - exon.getStart() + 1;
					if (pos < l) {
						int left_length = length - l + pos;
						the_read.getPositions().add(exon.getStart() + pos);
						while (left_length > 0) {
							j++;
							j %= script.getExons().size();
							if (exon.getEnd() + 1 == script.getExon(j).getStart()) {
								exon = script.getExon(j);
							}
							else {
								the_read.getPositions().add(exon.getEnd());
								exon = script.getExon(j);
								the_read.getPositions().add(exon.getStart());
							}
							left_length -= exon.getEnd() - exon.getStart() + 1;
						}
						the_read.getPositions().add(left_length + exon.getEnd());
						break;
					}
					else {
						pos -= l;
					}
				}
				the_read.setId(script.getId());
				the_read.setCirc_flag(true);
				the_read.setPeak_flag(true);
				out.add(the_read);
			}
			for (int i = 0; i < count_ip; i++) {
				int pos = Method.randIntReach(start_index, end_index);
				if (pos < 0) {
					pos += length;
				}
				Read the_read = new Read();
				for (int j = 0; ; j++) {
					j %= script.getExons().size();
					Exon exon = script.getExon(j);
					int l = exon.getEnd() - exon.getStart() + 1;
					if (pos < l) {
						int left_length = length - l + pos;
						the_read.getPositions().add(exon.getStart() + pos);
						while (left_length > 0) {
							j++;
							j %= script.getExons().size();
							if (exon.getEnd() + 1 == script.getExon(j).getStart()) {
								exon = script.getExon(j);
							}
							else {
								the_read.getPositions().add(exon.getEnd());
								exon = script.getExon(j);
								the_read.getPositions().add(exon.getStart());
							}
							left_length -= exon.getEnd() - exon.getStart() + 1;
						}
						the_read.getPositions().add(left_length + exon.getEnd());
						break;
					}
					else {
						pos -= l;
					}
				}
				the_read.setId(script.getId());
				the_read.setCirc_flag(true);
				the_read.setPeak_flag(true);
				out_ip.add(the_read);
			}
		}
		else {
			end_index = script.getBase_sum() - fragment_length * 2 - 1;
			start_index = Method.randIntReach(start_index, end_index);
			end_index = start_index + fragment_length;
			for (int i = 0; i < count; i++) {
				int pos = Method.randIntReach(start_index, end_index);
				Read the_read = new Read();
				for (int j = 0; j < script.getExons().size(); j++) {
					Exon exon = script.getExon(j);
					int l = exon.getEnd() - exon.getStart() + 1;
					if (pos < l) {
						int left_length = length - l + pos;
						the_read.getPositions().add(exon.getStart() + pos);
						while (left_length > 0) {
							the_read.getPositions().add(exon.getEnd());
							j++;
							exon = script.getExon(j);
							the_read.getPositions().add(exon.getStart());
							left_length -= exon.getEnd() - exon.getStart() + 1;
						}
						the_read.getPositions().add(left_length + exon.getEnd());
						break;
					}
					else {
						pos -= l;
					}
				}
				the_read.setId(script.getId());
				the_read.setPeak_flag(true);
				out.add(the_read);
			}
			for (int i = 0; i < count_ip; i++) {
				int pos = Method.randIntReach(start_index, end_index);
				Read the_read = new Read();
				for (int j = 0; j < script.getExons().size(); j++) {
					Exon exon = script.getExon(j);
					int l = exon.getEnd() - exon.getStart() + 1;
					if (pos < l) {
						int left_length = length - l + pos;
						the_read.getPositions().add(exon.getStart() + pos);
						while (left_length > 0) {
							the_read.getPositions().add(exon.getEnd());
							j++;
							exon = script.getExon(j);
							the_read.getPositions().add(exon.getStart());
							left_length -= exon.getEnd() - exon.getStart() + 1;
						}
						the_read.getPositions().add(left_length + exon.getEnd());
						break;
					}
					else {
						pos -= l;
					}
				}
				the_read.setId(script.getId());
				the_read.setPeak_flag(true);
				out_ip.add(the_read);
			}
		}
		
		ArrayList<Integer> pos = new ArrayList<>();
		for (int i = 0; ; i++) {
			i %= script.getExons().size();
			Exon exon = script.getExon(i);
			int l = exon.getEnd() - exon.getStart() + 1;
			if (start_index < l) {
				int left_length = length * 2 - l + start_index;
				pos.add(exon.getStart() + start_index);
				while (left_length > 0) {
					i++;
					i %= script.getExons().size();
					if (exon.getEnd() + 1 == script.getExon(i).getStart()) {
						exon = script.getExon(i);
					}
					else {
						pos.add(exon.getEnd());
						exon = script.getExon(i);
						pos.add(exon.getStart());
					}
					left_length -= exon.getEnd() - exon.getStart() + 1;
				}
				pos.add(left_length + exon.getEnd());
				break;
			}
			else {
				start_index -= l;
			}
		}
		Collections.sort(pos);
		Bed12 report = new Bed12();
		report.setBlock_count(pos.size() / 2);
		report.setChr(script.getExon(0).getChr_symbol());
		report.setEnd(pos.get(pos.size() - 1));
		report.setName(script.getId());
		report.setStart(pos.get(0));
		report.setStrand(script.getStrand());
		report.setThick_end(report.getEnd());
		report.setThick_start(report.getStart());
		for (int i = 0; i < pos.size(); i += 2) {
			report.getBlock_starts().add(pos.get(i) - report.getStart());
			report.getBlock_sizes().add(pos.get(i + 1) - pos.get(i));
		}
		return report;
	}
	
	private void getScriptReads(ArrayList<String> out, Transcript script, int count, StringBuffer id, int id_fix, StringBuffer base_seq, int start_index, int end_index, int length) {
		HashMap<Integer, Integer> count_id = new HashMap<Integer, Integer>();
		while (count > 0) {
			id.setLength(id_fix);
			base_seq.setLength(0);
			int rand = Method.randIntReach(start_index, end_index);
			if (count_id.containsKey(rand)) {
				count_id.put(rand, count_id.get(rand) + 1);
			}
			else {
				count_id.put(rand, 1);
			}
			id.append(count_id.get(rand));
			id.append(':');
			for (int i=0; i < script.getExons().size(); i++) {
				Exon exon = script.getExon(i);
				int l = exon.getEnd() - exon.getStart() + 1;
				if (rand < l) {
					int left_length = this.buildFirstSplice(exon, rand, length, base_seq, id);
					int j = i;
					while (left_length > 0) {
						j++;
						left_length = this.addBoundSplice(base_seq, id, left_length, script.getExon(j));
					}
					break;
				}
				else {
					rand -= l;
				}
			}			
			out.add(id.toString());
			out.add(base_seq.toString());
			count--;
		}
	}
	
//	int buildFirstSplice(Transcript script, int offset, int length, StringBuffer base_seq, StringBuffer id) {
//		int left_length = 0;
//		if (length > 0) {
//			String exon_base_seq = exon.getBase_seq();
//			left_length = length + offset - exon_base_seq.length() ;
//			if (left_length > 0) {
//				base_seq.append(exon_base_seq.substring(offset));
//				id.append(exon.getChr_symbol());
//				id.append(':');
//				id.append(exon.getStart() + offset);
//				id.append('-');
//				id.append(exon.getEnd());
//			}
//			else {
//				base_seq.append(exon_base_seq.substring(offset, offset + length));
//				id.append(exon.getChr_symbol());
//				id.append(':');
//				id.append(exon.getStart() + offset);
//				id.append('-');
//				id.append(exon.getStart() + offset + length - 1);
//			}
//		}
//		return left_length;
//	}
	
	int buildFirstSplice(Exon exon, int offset, int length, StringBuffer base_seq, StringBuffer id) {
		int left_length = 0;
		if (length > 0) {
			String exon_base_seq = exon.getBase_seq();
			left_length = length + offset - exon_base_seq.length() ;
			if (left_length > 0) {
				base_seq.append(exon_base_seq.substring(offset));
				id.append(exon.getChr_symbol());
				id.append(':');
				id.append(exon.getStart() + offset);
				id.append('-');
				id.append(exon.getEnd());
			}
			else {
				base_seq.append(exon_base_seq.substring(offset, offset + length));
				id.append(exon.getChr_symbol());
				id.append(':');
				id.append(exon.getStart() + offset);
				id.append('-');
				id.append(exon.getStart() + offset + length - 1);
			}
		}
		return left_length;
	}
	
	int addBoundSplice(StringBuffer base_seq, StringBuffer id, int length, Exon next_exon) {
		int left_length = 0;
		if (length > 0) {
			left_length = length - next_exon.getBase_seq().length();
			if (left_length > 0) {
				base_seq.append(next_exon.getBase_seq());
				appendID(id, next_exon.getStart(), next_exon.getEnd());
			}
			else {
				base_seq.append(next_exon.getBase_seq().substring(0, length));
				appendID(id, next_exon.getStart(), next_exon.getStart() + length - 1);
			}
		}
		return left_length;
	}
	
	void addAllBoundForwardSplice(StringBuffer base_seq, StringBuffer id, int length, int exon_index, ArrayList<Exon> exon_list){
		if (length > 0) {
			int next_index = exon_list.size() - 1;
			if (exon_index < next_index) {
				next_index = Method.randIntReach(exon_index + 1, next_index);
			}
			else {
				System.out.println("Warning: Last exon forward splice");
				id.append('L');
			}
			Exon next_exon = exon_list.get(next_index);
			int left_length = length - next_exon.getBase_seq().length();
			if (left_length > 0) {
				base_seq.append(next_exon.getBase_seq());
				appendID(id, next_exon.getStart(), next_exon.getEnd());
				addAllBoundForwardSplice(base_seq, id, left_length, next_index, exon_list);
			}
			else {
				base_seq.append(next_exon.getBase_seq().substring(0, length));
				appendID(id, next_exon.getStart(), next_exon.getStart() + length - 1);
			}
		}
	}
	
	private void appendID(StringBuffer id, int start, int end) {
		id.append('_');
		id.append(start);
		id.append('-');
		id.append(end);
	}
}
