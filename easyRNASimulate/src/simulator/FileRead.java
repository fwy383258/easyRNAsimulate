package simulator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class FileRead {
	
	public BufferedReader openFile(String file_name) {
		BufferedReader reader=null;
		File read_file=new File(file_name);
		try {
			reader=new BufferedReader(new FileReader(read_file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		return reader;
	}
	
	public void closeFile(BufferedReader reader) {
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
