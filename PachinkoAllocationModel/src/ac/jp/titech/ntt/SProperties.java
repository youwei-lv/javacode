package ac.jp.titech.ntt;

import java.io.IOException;
import java.io.StringReader;
import java.util.Properties;

/**
 * SProperties is derived from java.io.Properties, with the extension that it
 * defines the method of loading properties from a String array, which might
 * come from the arguments passed to the main function.
 * 
 * @author Youwei Lu
 */
public class SProperties extends Properties {

	private static final long serialVersionUID = 6957966910702529508L;

	/**
	 * Loads properties from a String array, which can be the arguments passed
	 * to the main function, as an additional member of the methods of loading
	 * properties from different sources defined in the superclass.
	 * <p>
	 * 
	 * Despite the main purpose to conveniently read the passed arguments of the
	 * main function, it can be generally used whenever the properties are
	 * stored in a String array, yet with the restriction that the resulted
	 * string obtained from concatenating the array strings with '\n' could be
	 * interpreted by the load(java.io.Reader) method in the superclass.
	 * <p>
	 * 
	 * "-" and "--" are often used to indicate the option names in command
	 * lines; however, they will be automatically stripped so that only the
	 * option names are stored. If other indicators are employed, refer to
	 * {@link #loadFromStringArray(String[], String)}.
	 * 
	 * @param args
	 *            String array storing the properties to be read
	 */
	public void load(String[] args) {
		loadFromStringArray(args, "(^-)|(\\s+-)|(^--)|(\\s+--)");
	}

	/**
	 * Works in the same way as {@link #load(String[])}, except that the
	 * indicators for the option names can be specified.
	 * 
	 * @param args
	 *            String array storing the properties to be read
	 * @param optionNameIndicator
	 *            Regular expression to identify the option names
	 */
	public void load(String[] args, String optionNameIndicator) {
		loadFromStringArray(args, optionNameIndicator);
	}

	public void checkNullKeys(String[] keys) {
		for(String k: keys) {
			if(!this.containsKey(k)) {
				throw new RuntimeException(k+" value is missing.");
			}
		}
	}
	
	private void loadFromStringArray(String[] args, String optionNameIndicatorPattern) {
		StringBuilder sb = new StringBuilder();
		for (String arg : args) {
			sb.append(arg).append("\n");
		}
		String commandlineProperties = sb.toString().replaceAll(optionNameIndicatorPattern, "\n");
		StringReader read = new StringReader(commandlineProperties);
		try {
			this.load(read);
		} catch (IOException e) {
			e.printStackTrace();
		}
		read.close();
	}

	static public void main(String[] args) {
		SProperties props = new SProperties();
		props.load(args);
		props.list(System.out);
		System.out.flush();
	}

}
