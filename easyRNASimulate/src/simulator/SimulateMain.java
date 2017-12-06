package simulator;


public class SimulateMain {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ReadInputs read_args = new ReadInputs(args);
		int c = read_args.completeArgs();
		if (c != 0) {
			if (c < read_args.getFULL_NOTE()) {
				System.out.println("Lack parameter: ");
				read_args.printNote(c);
			}
		}
		else {
			SimulateModes.run(read_args);
		}
		
		System.out.println("END");
	}
	


}
