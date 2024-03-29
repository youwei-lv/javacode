package eng;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Properties;

import edu.stanford.nlp.ling.CoreAnnotations.LemmaAnnotation;
import edu.stanford.nlp.ling.CoreAnnotations.SentencesAnnotation;
import edu.stanford.nlp.ling.CoreAnnotations.TextAnnotation;
import edu.stanford.nlp.ling.CoreAnnotations.TokensAnnotation;
import edu.stanford.nlp.ling.CoreLabel;
import edu.stanford.nlp.pipeline.Annotation;
import edu.stanford.nlp.pipeline.StanfordCoreNLP;
import edu.stanford.nlp.util.CoreMap;
import util.SProperties;

/**
 * An object of BaseFormConverter convert English texts to texts of words in the base form by using
 * <a href="http://nlp.stanford.edu/software/corenlp.shtml">Standford Core NLP Tools</a>.
 * <p>
 * The input can be a single file or a directory. When the input is a single file, the converting
 * results will be written to the specified output file. Note the specified output file must not be 
 * an existing directory; otherwise, an error will occur making the processing fail.<br>
 * If the input is an existing directory, the files in the directory will be processed individually
 * with the words in the base form written to a file in the specified output directory with the same 
 * file name as the source file. It will also cause a failure if the specified output is 
 * an existing file.
 * <p>
 * 
 * @author Youwei Lu
 */
public class BaseFormConverter {
	/**
	 * Input file name which can be a single text file or a directory.
	 */
	String inputf;
	/**
	 * Output file name. If the input is a single file, it is the output file name.
	 * If the input is a directory, it becomes the directory name carrying the output files from
	 * processing files in the input directory individually.
	 */
	String outputf;
	/**
	 * Error file name, "error.log" by default.
	 */
	String errorlogfile;
	
	private boolean isdir;
	
	// a StanfordCoreNLP object, with POS tagging, lemmatization, NER, parsing, and coreference resolution 
	private StanfordCoreNLP pipeline;
	
	public BaseFormConverter(String inputf, String outputf, String errorlogfile) {
		this.inputf = inputf;
		this.outputf = outputf;
		this.errorlogfile = errorlogfile;
		File output = new File(outputf);
		if(new File(inputf).isDirectory()){
			this.isdir = true;
			if(output.exists() && output.isFile()) {
				throw new RuntimeException("The output directory cannot be created");
			}
			else if(!output.exists()) {
				output.mkdir();
			}
		}
		else {
			if(output.exists() && output.isDirectory()) {
				throw new RuntimeException("Directory "+outputf+" is already used.");
			}
			this.isdir = false;
		}
		Properties props = new Properties();
		props.setProperty("annotators", "tokenize, ssplit, pos, lemma");
		//props.setProperty("annotators", "tokenize, ssplit, pos, lemma, ner, parse, dcoref");
		pipeline = new StanfordCoreNLP(props);
	}
	
	public BaseFormConverter(String inputf, String outputf) {
		this(inputf, outputf, "error.log");
	}
	
	public void doconvert() {
		PrintWriter err;
		File outputfile = null;
		BufferedWriter writer = null;
		BufferedReader reader = null;
		String line = null;
		List<CoreMap> sentences = null;
		try {
			err = new PrintWriter(this.errorlogfile);
		} catch (FileNotFoundException e1) {
			throw new RuntimeException("Cannot create log file.");
		}
		if(isdir) {
			File inputdir = new File(this.inputf);
			File[] directoryListing = inputdir.listFiles();
			for (File child : directoryListing) {
				System.out.println("processing file "+child.getName());
				try {
					outputfile = new File(outputf, child.getName());
					writer = new BufferedWriter(new FileWriter(outputfile));
					reader = new BufferedReader(new FileReader(child));
					while((line=reader.readLine()) != null) {
						sentences = convertStringToSentences(line);
						outputSentencesInBasicForm(sentences, writer);
						writer.newLine();
					}
				} catch (IOException e) {
					e.printStackTrace(err);
				} finally {
					try {
						if(writer != null) writer.close();
						if(writer != null) reader.close();
					} catch(IOException e) {
						e.printStackTrace(err);
					}
				}
			}
		}
		else {
			try {
				reader = new BufferedReader(new FileReader(this.inputf));
				writer = new BufferedWriter(new FileWriter(this.outputf));
				while((line=reader.readLine()) != null) {
					sentences = convertStringToSentences(line);
					outputSentencesInBasicForm(sentences, writer);
					writer.newLine();
				}
			} catch (IOException e) {
				e.printStackTrace(err);
			} finally {
				try {
					if(writer != null) writer.close();
					if(writer != null) reader.close();
				} catch(IOException e) {
					e.printStackTrace(err);
				}
			}
		}
	}
	
	static void outputSentencesInBasicForm(List<CoreMap> sentences, BufferedWriter out) throws IOException {
		@SuppressWarnings("unused")
		String word, lemma;
		for(CoreMap sentence: sentences) {
			for (CoreLabel token: sentence.get(TokensAnnotation.class)) {
		        // this is the text of the token
		        word = token.get(TextAnnotation.class);
		        // this is the lemma of the token
		        lemma = token.get(LemmaAnnotation.class);
		        out.write(lemma+" ");
		    }
		}
		out.flush();
	}
	
	public List<CoreMap> convertStringToSentences(String text) {
	    // create an empty Annotation just with the given text
	    Annotation document = new Annotation(text);
	    
	    // run all Annotators on this text
	    pipeline.annotate(document);
	    
	    // these are all the sentences in this document
	    // a CoreMap is essentially a Map that uses class objects as keys and has values with custom types
	    List<CoreMap> sentences = document.get(SentencesAnnotation.class);
	    
	    return sentences;
	}

	public static void main(String[] args) {
		SProperties props = new SProperties();
		props.load(args);
		
		String input, output;
		if((input=props.getProperty("input"))==null) {
			throw new RuntimeException("Use -input=<filename> option to specify the input file/directory.");
		}
		if((output=props.getProperty("output"))==null) {
			throw new RuntimeException("Use -output=<filename> option to specify the output file.");
		}
		
		BaseFormConverter converter = new BaseFormConverter(input, output);
		converter.doconvert();
	}

}
