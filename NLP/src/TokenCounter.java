import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * A TokenCounter object counts the appearance times of different tokens in
 * files.
 * 
 * @author Youwei Lu
 */
public class TokenCounter {

	/**
	 * Specified counting target, which may be a file name or a directory. When
	 * it is a single file, the tokens in the file will be counted; however,
	 * when is a directory, tokens in all files in the directory will be
	 * counted.
	 */
	String inputf;
	/**
	 * The name of file to which the counting result will be written.
	 */
	String countf;
	/**
	 * A regular expression for matched tokens to be filtered out from counting
	 */
	String ignorePattern;
	/**
	 * A counting map between tokens and their frequencies
	 */
	private Map<String, Integer> count;
	/**
	 * Constant standing for sorting in ascending order.
	 */
	static final int ASC = 0;
	/**
	 * Constant standing for sorting in descending order.
	 */
	static final int DESC = 1;

	public TokenCounter(String inputf, String countf) {
		this.inputf = inputf;
		this.countf = countf;
		count = new HashMap<String, Integer>();
	}

	/**
	 * Increase the frequency of a specific token by one.
	 */
	public void addToken(String token) {
		if (count.containsKey(token)) {
			count.put(token, count.get(token) + 1);
		} else {
			count.put(token, 1);
		}
	}

	/**
	 * Output a list of token-frequency pairs which are sorted according to
	 * their frequencies in the ascending or descending order.
	 * 
	 * @param order Mode of sorting, which is either {@link #ASC} or {@link #DESC}.
	 * @return A sorted list of token-frequency pairs.
	 */
	public List<Map.Entry<String, Integer>> sortByFrequency(int order) {
		List<Map.Entry<String, Integer>> list = new LinkedList<>(count.entrySet());
		if (order == ASC) {
			Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {
				@Override
				public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
					return o1.getValue().compareTo(o2.getValue());
				}
			});
		} else if (order == DESC) {
			Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {
				@Override
				public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
					return -o1.getValue().compareTo(o2.getValue());
				}
			});
		} else {
			throw new RuntimeException("Unknown sorting order.");
		}
		return list;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}