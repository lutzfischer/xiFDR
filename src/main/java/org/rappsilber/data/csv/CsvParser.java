/*
 * Copyright 2015 Lutz Fischer <lfischer at staffmail.ed.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.rappsilber.data.csv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.rappsilber.utils.UpdatableChar;
import org.rappsilber.utils.UpdateableInteger;

/**
 * implements a generic csv-parser, that enables a to access fields by name - including alternatives to these names.
 * @author lfischer
 */
public class CsvParser {
    /** array, that contains the header of the file **/
    private String[]                    m_header;
    /** list of list, containing alternatives for column-names **/
    private ArrayList<HashSet<String>>  m_headerAlternatives = new ArrayList<HashSet<String>>();
    /** maps the column-name to the column-index **/
    private HashMap<String,Integer>     m_headerToColumn = new HashMap<String, Integer>();
    /** the field delimeter **/
    private String                      m_delimiter = ",";
    /** quote character **/
    private String                      m_quote = "\"";
    /** just a short-cut to the duplicate quotes - for having the quote character within quotes **/
    private String                      m_doublequote = m_quote+m_quote;
    
    /** a pattern, that matches the values in a cell. it provides two groups - first for quoted cell-values and second for unquoted cell-values **/
    private Pattern                     m_cellValue;
                                        //Pattern.compile(m_delimiter +"?\\s*(?:\"((?:[^\"]*(?:\"\")?)*)\"|([^"+ m_delimiter +"]*))\\s*"+ m_delimiter + "?");
    /** the maximum number of columns found sofar in the csv-file **/
    private int                         m_maxColumns=0;
    /** how many columns where in the header-row **/
    private int                         m_countHeader;
    /** the current line, that is parsed **/
    private String                      m_currentLine;
    /** the current line split up into values **/
    private String[]                    m_currentValues;
    /** how many values where found in the current line **/
    private int                         m_foundValues=0;
    /** the number of the current line **/
    private int                         m_lineNumber=0;
    /** the reader for the file - if it reads from a file*/
    private BufferedReader              m_input;
    /** the file that is read - if a file is read **/
    private File                        m_inputFile;
    /** if we should guess the delimiter - what are the candidates**/
    protected static char[]             TEST_DELIMITERS = {',','\t','|'};
    /** if we should guess the quote chars - what are the candidates**/
    protected static char[]             TEST_QUOTES = {'\'','"'};
    /** a pattern defining the values, that would interpreted as true */
    public final static Pattern      ISTRUE = Pattern.compile("^(T|1|-1|TRUE|Y|YES|\\+)$",Pattern.CASE_INSENSITIVE);
    /** a pattern defining the values, that would interpreted as false */
    public final static Pattern      ISFALSE = Pattern.compile("^(F|0|FALSE|N|NO|\\-)?$",Pattern.CASE_INSENSITIVE);
    
    
    
    /** 
     * if a field is missing this is the value it gets 
     * As it is an independent object we can just determine if a field where missing by comparing (=)
     */
    protected String                      MISSING_FIELD = new String("");

    /**
     * default constructor
     */
    public CsvParser() {
        setPattern(m_delimiter, m_quote);
    }

    /** 
     * constructor defining the field delimiter
     * @param delimiter 
     */
    public CsvParser(char delimiter) {
        setPattern(Character.toString(delimiter), m_quote);
    }
    
    /** 
     * constructor defining quote character
     * @param delimiter 
     * @param quote 
     */
    public CsvParser(char delimiter, char quote) {
        setDelimiter(delimiter);
        setQuote(quote);
    }
    
    /**
     * create the pattern used to split up rows according to the delimeter and quote characters
     * @param d
     * @param q 
     */
    protected void setPattern(String d, String q) {
        m_delimiter = d;
        m_quote = q;
        m_doublequote = m_quote + m_quote;
        if (d.matches("\\t")) {
            m_cellValue = Pattern.compile(" *(?:"+q+"((?:[^"+q+"]*(?:"+q+q+")?)*)"+q+"|([^"+ d +"]*[^ "+ d +"])|) *"+ d + "?");
        } else if (d.matches(" ")) {
            m_cellValue = Pattern.compile("\\t*(?:"+q+"((?:[^"+q+"]*(?:"+q+q+")?)*)"+q+"|([^"+ d +"]*[^\\t"+ d +"])|)\\t*"+ d + "?");
        } else
            m_cellValue = Pattern.compile("\\s*(?:"+q+"((?:[^"+q+"]*(?:"+q+q+")?)*)"+q+"|([^"+ d +"]*[^\\s"+ d +"])|)\\s*"+ d + "?");
    }

    
    /**
     * create the pattern used to split up rows according to the delimeter and quote characters
     * @param d
     * @param q 
     */
    public void setPattern(char d, char q) {
        setPattern(Character.toString(d), Character.toString(q));
    }

    
    /**
     * uses the given pattern to split rows into fields
     * 
     * @param p pattern used to separate the rows into individual fields 
     */
    protected void setPattern(Pattern p) {
        m_cellValue = p;
        m_delimiter = "";
        m_quote = "";
        m_doublequote = "";
    }

    /**
     * set the character that separates two adjunct fields
     * @param d 
     */
    public void setDelimiter(char d) {
        String delim = Character.toString(d);
        if (!("|.\\()[]^$".contains(delim)))
            setPattern(delim, m_quote);
        else
            setPattern("\\" + delim, m_quote);
    }

    /**
     * set the character that can be used to surround a field-value if it contains the delimiter-character
     * @param q 
     */
    public void setQuote(char q) {
        setPattern(m_delimiter, Character.toString(q));
    }

    /**
     * return the delimiter-character
     * @return 
     */
    public String getDelimiter() {
        return m_delimiter;
    }

    /**
     * return the quote-character
     * @return 
     */
    public String getQuote() {
        return m_quote;
    }

    /**
     * split a line according to the current configuration
     * @param line
     * @return 
     */
    public ArrayList<String> splitLine(String line) {
        return splitLine(line, new UpdateableInteger(0));
    }
    
    /**
     * split a line according to the current configuration
     * @param line
     * @param countQuoted returns the number of quoted fields
     * @return 
     */
    public ArrayList<String> splitLine(String line, UpdateableInteger countQuoted) {
        ArrayList<String> ret = new ArrayList<String>(m_maxColumns);
        Matcher m = m_cellValue.matcher(line);
        
        int c = 0;
        while (m.find()) {

            if (m.group(1) != null) {
                ret.add(m.group(1).replace(m_doublequote, m_quote));
            } else {
                ret.add(m.group(2));
                countQuoted.value++;
            }
            c++;
        }
//        if (c>m_maxColumns)
//            m_maxColumns = c;
        ret.remove(ret.size() - 1);
        return ret;
    }
    
    /**
     * defines the csv-header according to the given line
     * @param line 
     */
    public void setHeader(String line) {
        ArrayList<String> hc = splitLine(line);
        // first set the header array
        String[] header = new String[hc.size()];
        hc.toArray(header);
        setHeader(header);
    }

    /**
     * Sets the column-headers to the given names
     * @param header
     */
    public void setHeader(String[] header) {
        m_header = header;
        if (header == null) {
            getHeaderToColumn().clear();
        } else {
            for (int i = 0; i< m_header.length; i++) {
                String h = m_header[i];
                if (h==null) 
                    h="Column "+Integer.toString(i);
                else 
                    h = h.trim();
                m_header[i]=h;
                getHeaderToColumn().put(h , i);
                String hl = h.toLowerCase();
                if (!h.contentEquals(hl) && getHeaderToColumn().get(h.toLowerCase()) == null)
                    getHeaderToColumn().put(hl , i);
                // now we look, if we have some alternative names for the column
                for (HashSet<String> names : getHeaderAlternatives()) {
                    if (names.contains(h) || names.contains(hl)) {
                        for (String ha : names) {
                            getHeaderToColumn().put(ha , i);
                        }
                        break;
                    }
                }
            }
        }
    }
    
    /**
     * removes any column-name
     */
    public void cleanHeader() {
        m_header = null;
        getHeaderToColumn().clear();
    }

    
    /**
     * returns an array with names of each column
     * @return 
     */
    public String[] getHeader() {
        return m_header;
    }
    
    
    /**
     * returns the name of the given column
     * @param column
     * @return 
     */
    public String getHeader(int column) {
        if (m_header == null || column >= m_header.length) {
            return Integer.toString(column);
        }
        String ret = m_header[column];
        if (ret == null || ret.trim().length() == 0)
            return Integer.toString(column);
        return m_header[column];
    }
    
    /**
     * declare an alternative for the given column name. 
     * If the file contains a column with either of these names - both names can 
     * be used to refer to the column.
     * @param name
     * @param alternative 
     */
    public void setAlternative(String name, String alternative) {
        boolean isSet=false;
        for (HashSet<String> names : getHeaderAlternatives()) {
            if (names.contains(name)) {
                names.add(alternative);
                Integer col = getColumn(name);
                if (col == null) {
                    col = getColumn(alternative);
                    if (col != null)
                        for (String n : names) {
                            getHeaderToColumn().put(n, col);
                        }
                } else 
                    getHeaderToColumn().put(alternative, col);
                
                isSet = true;
                break;
            }
            if (names.contains(alternative)) {
                String reference = names.iterator().next();
                names.add(name);
                Integer col = getColumn(reference);
                if (col != null)
                    getHeaderToColumn().put(name, col);
                
                isSet = true;
                break;
            }
        }
        if (!isSet) {
            HashSet<String> names = new HashSet<String>(2);
            names.add(name);
            names.add(alternative);
            getHeaderAlternatives().add(names);
            if (getColumn(name) != null)
                getHeaderToColumn().put(alternative, getColumn(name));
            else if (getColumn(alternative) != null)
                getHeaderToColumn().put(name, getColumn(alternative));
        }
    }


    /**
     * sets the next line as header
     * @throws IOException 
     */
    public void readHeader() throws IOException {
        next();
        setHeader(getCurrentLine());
    }
    
    /**
     * define the String as the current line.
     * @param line 
     */
    public void setCurrentLine(String line) {
        m_lineNumber++;
        m_currentLine = line;
        m_currentValues = splitLine(line).toArray(new String[0]);
        for (int c = 0; c<m_currentValues.length;c++)
            if (m_currentValues[c] == null)
                m_currentValues[c] = MISSING_FIELD;
        int l = m_currentValues.length;
        if (l > m_maxColumns || m_maxColumns == 0) {
            if (m_maxColumns != 0) {
                Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "line " + 
                                            m_lineNumber + 
                                            " contains more fields than ever observed!\n" +
                                            line);
            }
            
            m_maxColumns = l;
            if (m_header != null) {
                String[] dheader = java.util.Arrays.copyOf(m_header, m_maxColumns);
                java.util.Arrays.fill(dheader, m_header.length, m_maxColumns, MISSING_FIELD);
                m_header = dheader;
            }
            
//            m_currentValues = new String[m_maxColumns];
            if (m_header != null) {
                if (l<m_header.length) {
                    Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "line " + m_lineNumber + " less fields than the header\n" +
                                                                    "Assume the missing values are at the end and will be empty");
                }
            }
        }
        
        m_foundValues = l;
    }

    /**
     * opens the file assuming it contains no column-headers
     * @param f
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void openFile(File f) throws FileNotFoundException, IOException {
        openFile(f, false);
    }
    
    /**
     * opens the file
     * @param f the file to be opened
     * @param hasHeader should the first row be interpreted as header
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void openFile(File f, boolean hasHeader) throws FileNotFoundException, IOException {
        m_inputFile = f;
        m_input = new BufferedReader(new FileReader(f));
        if (hasHeader)
            readHeader();
    }
    
    public void openURI(URI f, boolean hasHeader) throws FileNotFoundException, IOException {
        if (f.getScheme().equals("file")) {
            openFile(new File(f), hasHeader);
        } else {
            URL url = new URL(f.toString());
            URLConnection connection = url.openConnection();
            InputStream is = connection.getInputStream();    
            m_inputFile = null;
            m_input = new BufferedReader(new InputStreamReader(is));
            if (hasHeader)
                readHeader();
        }
    }
    
    public void openURI(URI f) throws FileNotFoundException, IOException {
        openURI(f, false);
    }
    
    /**
     * takes the current line, and sets it as header
     */    
    public void setCurrentLineAsHeader() {
        setHeader(getCurrentLine());
    }
    
    /**
     * creates a new csv-parser for the given file
     * @param f the file to parse
     * @param delimiter field-delimiter
     * @param quote quote character
     * @return a csv-parser for the given file
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static CsvParser getCSV(File f, char delimiter, char quote) throws FileNotFoundException, IOException {
        CsvParser csv = new CsvParser(delimiter, quote);
        csv.openFile(f);
        return csv;
    }

    /**
     * Opens a csv-parser for the given file and tries to predict the quote 
     * character and the field-delimiter based on the first number of lines.<br/>
     * After it reads the first lines it will reopen the file with a configured 
     * csv-parser.<br/>
     * <b><i>This means, that the file must be re-readable e.g no pipe</i></b>
     * 
     * @param f the file to open
     * @param guessLines how many lines to read in for guessing the parameters
     * @return a {@link CavParser} for the file
     * @throws FileNotFoundException
     * @throws IOException
     */
    public static CsvParser guessCsv(File f, int guessLines) throws FileNotFoundException, IOException {
        return guessCsv(f, guessLines, TEST_DELIMITERS, TEST_QUOTES);
    }
    
    /**
     * Opens a csv-parser for the given file and tries to predict the quote 
     * character and the field-delimiter based on the first number of lines.<br/>
     * After it reads the first lines it will reopen the file with a configured 
     * csv-parser.<br/>
     * <b><i>This means, that the file must be re-readable e.g no pipe</i></b>
     * 
     * @param f the file to open
     * @param guessLines how many lines to read in for guessing the parameters
     * @param delimiters an array with chars that should be considered as delimiters
     * @param quotes an array with chars that should be considered as quote-characters
     * @return a {@link CavParser} for the file
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static CsvParser guessCsv(File f, int guessLines, char[] delimiters, char[] quotes) throws FileNotFoundException, IOException {
      
        UpdatableChar delimiter = new UpdatableChar();
        UpdatableChar quote = new UpdatableChar();
        Boolean unique = guessDelimQuote(f, guessLines, delimiters, quotes, delimiter, quote);
                
        
        if (unique == null) {
            return null;
        }
        
        CsvParser csv = new CsvParser(delimiter.value, quote.value);
        csv.openFile(f);
        return csv;
    }
    
    /**
     * opens the file and reads a guesslines number of rows and tries to guess what are the correct delimiter and quote characters
     * @param f
     * @param guessLines
     * @param delimiter returns the guest delimiter char
     * @param quote returns the guest quote char
     * @return true if there was one setting of delimiter and quote that was better than any other
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static Boolean guessDelimQuote(File f, int guessLines, UpdatableChar delimiter, UpdatableChar quote) throws FileNotFoundException, IOException {
        return guessDelimQuote(f, guessLines, TEST_DELIMITERS, TEST_QUOTES, delimiter, quote);
    }
    
    /**
     * opens the file and reads a guesslines number of rows and tries to guess 
     * what are the correct delimiter and quote characters from the provided 
     * choices
     * @param f
     * @param guessLines
     * @param delimiters set of delimiter characters from which the correct one 
     * should be guessed.
     * @param quotes set of quote characters from which the correct one 
     * should be guessed.
     * @param delimiter returns the guest delimiter char
     * @param quote returns the guest quote char
     * @return true if there was one setting of delimiter and quote that was better than any other
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static Boolean guessDelimQuote(File f, int guessLines, char[] delimiters, char[] quotes, UpdatableChar delimiter, UpdatableChar quote) throws FileNotFoundException, IOException {
        ArrayList<String> allTestLines = new ArrayList<String>(guessLines);
        int lc = 0;
        String line;
        BufferedReader br = new BufferedReader(new FileReader(f));
        
        while ((line = br.readLine()) != null && lc++<guessLines) {
            allTestLines.add(line);
        }
        
        if (allTestLines.isEmpty()) {
            delimiter.value = delimiters[0];
            quote.value = quotes[0];
            return false;
        }
        
        
        int[] values = new int[allTestLines.size()];
        char delimChar = ' ';
        char quoteChar = ' ';
        int minFieldCount = 0;
        int maxFieldCount = 0;
        boolean unique = false;
        int quoteFound = 0;
        
        dloop: for (char d : delimiters) {
            int dMinFieldCount = 0;
            int dMaxFieldCount = 0;
            int dQuoteFound = 0;
            char dQuoteChar = ' ';
            
            qloop: for (char q : quotes) {
                int thisquoteFound = 0;
                CsvParser testCsv = new CsvParser(d, q);
                // how many fields would we have in the first line
                int c = testCsv.splitLine(allTestLines.get(0)).size();
                int qMinFieldCount = c;
                int qMaxFieldCount = c;
                // if it is only one field - we cant conclude anything -ignored
//                if (c  == 1 )
//                    continue qloop;
                
                UpdateableInteger countQuoted = new UpdateableInteger(0);
                for (String l : allTestLines) {
                    c=testCsv.splitLine(l,countQuoted).size();
                            
                    // if not all lines have the same number of fields we ignore this combination
                    if (qMinFieldCount > c)
                        qMinFieldCount = c;
                    else if (qMaxFieldCount < c)
                        qMaxFieldCount = c;
                    thisquoteFound += l.length() - l.replace(""+q, "").length();
                }
                
                if (dMinFieldCount == 0) {
                    dMinFieldCount = qMinFieldCount;
                    dMaxFieldCount = qMaxFieldCount;
                    dQuoteChar = q;
                    unique = true;
                    dQuoteFound = thisquoteFound;
                } else if (qMinFieldCount>=dMinFieldCount || 
                        (qMinFieldCount == dMinFieldCount && qMaxFieldCount < dMaxFieldCount) ||
                        (qMinFieldCount == dMinFieldCount && qMaxFieldCount == dMaxFieldCount && dQuoteFound<thisquoteFound)) {
                    dMinFieldCount = qMinFieldCount;
                    dMaxFieldCount = qMaxFieldCount;
                    dQuoteChar = q;
                    unique = true;
                    dQuoteFound = thisquoteFound;
                } else if (qMinFieldCount == dMinFieldCount && qMaxFieldCount == dMaxFieldCount && dQuoteFound == thisquoteFound) {
                    unique = false;
                }
            }
            if (dMinFieldCount>=minFieldCount || 
                (dMinFieldCount == minFieldCount && dMaxFieldCount > maxFieldCount)) {
                quoteChar  = dQuoteChar;
                delimChar = d;
                minFieldCount = dMinFieldCount;
                maxFieldCount = dMaxFieldCount;
            }
            
        }
        
        if (maxFieldCount == 0) {
            return null;
        }
        
        delimiter.value = delimChar;
        quote.value = quoteChar;
        return unique;
        
    }
    
    /**
     * reads the next row from the file
     * @return true if there is actually a next row - or false if we already reached the end of the file
     * @throws java.io.IOException
     */
    public boolean next() throws IOException {
        if (m_input == null)
            throw new IOException("No file is currently open for reading");
        String line = m_input.readLine();
        while (line != null && line.isEmpty())
            line = m_input.readLine();
        
        if (line  == null)
            return false;
        setCurrentLine(line);
        return true;
    }
    
    /**
     * does the current row only contains empty fields
     * @return 
     */
    public boolean onlyEmptyFields() {
        for (int i = 0; i< m_foundValues; i++) {
            if (!m_currentValues[i].isEmpty())
                return false;
        }
        return true;
    }
    
    /**
     * adds a new unnamed column to the csv-file
     * @return the index of the column
     */
    public int addColumn() {
        return m_maxColumns++;
    }
    
    /**
     * adds a new named column to the csv-file
     * @param name the name of the column
     * @return the index of the new column
     */
    public int addColumn(String name) {
        int col = addColumn();
        if (m_header != null) {
            String[] dheader = java.util.Arrays.copyOf(m_header, col+1);
            dheader[col] = name;
            m_header = dheader;
        }
        getHeaderToColumn().put(name, col);
        if (!name.contentEquals(name.toLowerCase()) && getHeaderToColumn().get(name) == null)
            getHeaderToColumn().put(name.toLowerCase(), col);
        return col;
    }

    
    /**
     * Tests if a column of the given name exists and if it does it returns the index of the existing column.
     * If no column of the given name exists a new column is added and the index of the new column is returned.
     * @param name the name of the column
     * @return the index of the existing or newly created column
     */
    public synchronized int addGetColumn(String name) {
        Integer c = getColumn(name);
        if (c == null) {
            int ret = addColumn(name);
            return ret;
        }
        return c;
    }

    /**
     * is it an empty row (no values and also no delimiters found)
     * @return 
     */
    public boolean onlyMissingFields() {
        return m_foundValues == 0;
    }
    
    /**
     * does the current row only contains empty fields
     * @return 
     */
    public boolean isEmpty() {
        return onlyEmptyFields();
    }

    /**
     * does the current row has an no value for the given field?
     * @param field name of the field/column
     * @return 
     */
    public boolean isEmpty(String field) {
        Integer column = getHeaderToColumn().get(field);
        if (column == null)
            return true;
        return m_currentValues[column].isEmpty();
    }

    /**
     * Does the current row has an no value for the given field?
     * <p>Also missing fields are considered empty</p>
     * @param field index of the field/column
     * @return true if the field is empty or missing; false otherwise
     */
    public boolean isEmpty(int field) {
        if (field > m_maxColumns)
            return true;
        return m_currentValues[field].isEmpty();
    }
    
    /**
     * Checks whether the given field is missing as opposed to be empty 
     * <br/> E.g. ,, would be both missing and empty but ,"", would be just empty.
     * @param field
     * @return 
     */
    public boolean isMissing(int field) {
        if (field > m_maxColumns)
            return true;
        return m_currentValues[field] == MISSING_FIELD;
    }
    
    /**
     * Checks whether the given field is missing as opposed to be empty 
     * <br/> E.g. ,, would be both missing and empty but ,"", would be just empty.
     * @param field
     * @return 
     */
    public boolean isMissing(String field) {
        Integer column = getHeaderToColumn().get(field);
        if (column == null)
            return true;
        return m_currentValues[column] == MISSING_FIELD;
    }
    
    /**
     * Returns the column id for the given name - or the ID of a column that has 
     * an accepted alternative name
     * @param name the name of the column
     * @return the id of the column - if matching column was found; null 
     * otherwise
     */
    public Integer getColumn(String name)  {
        Integer col = getHeaderToColumn().get(name);
        if (col == null && name.matches("[0-9]*")) {
            Integer name2col = Integer.parseInt(name);
            if (name2col < getMaxColumns())
                return name2col;
        }
        
        return col;
    }
    
    /**
     * returns the String-value of the given field
     * @param field
     * @return 
     */
    public String getValue(String field) {
        Integer column = getHeaderToColumn().get(field);
        if (column == null)
            return MISSING_FIELD;
        return m_currentValues[column];
    }
    /**
     * returns the String-value of the given field
     * @param field
     * @return 
     */
    public String getValue(int field) {
        if (field > m_maxColumns)
            return MISSING_FIELD;
        return m_currentValues[field];
    }
    /**
     * returns the String-value of the given field
     * @param field
     * @return 
     */
    public String getValue(Integer field) {
        if (field == null || field >= m_currentValues.length)
            return MISSING_FIELD;
        return m_currentValues[field];
    }
    
    /**
     * Returns the String-value of the given field
     * <p>If the field is missing then defaultValue is used to calculate a value
     * @param field
     * @param defaultValue
     * @return 
     */    
    public String getValue(Integer field, CSVValueCalc defaultValue) {
        if (field == null || field > m_maxColumns) {
            return defaultValue.getValue(this);
        }
        if (m_currentValues[field] == MISSING_FIELD) {
            String v = defaultValue.getValue(this);
            m_currentValues[field] = v;
            return v;
        }
        return m_currentValues[field];
    }
    
    /**
     * Returns the String-value of the given field
     * <p>If the field is missing then defaultValue is used to calculate a value
     * @param field
     * @param defaultValue
     * @return 
     */        
    public String getValue(String field, CSVValueCalc defaultValue) {
        Integer column = getHeaderToColumn().get(field);
        if (column == null)
            return defaultValue.getValue(this);
        
        return getValue(column, defaultValue);
    }    
    
    /**
     * Returns the value of the given field interpreted as Integer.
     * <p> Internally the String get converted via 
     * {@link Integer#parseInt(java.lang.String)}.
     * @param field
     * @return 
     */    
    public Integer getInteger(Integer field) {
        String v = getValue(field);
        if (v == MISSING_FIELD)
            return null;
        
        return Integer.parseInt(v);
    }
    
    /**
     * Returns the value of the given field interpreted as Integer.
     * <p> Internally the String get converted via 
     * {@link Integer#parseInt(java.lang.String)}.</p>
     * @param field
     * @param defaultValue  if the field is missing this is returned
     * @return 
     */    
    public Integer getInteger(Integer field, Integer defaultValue) {
        String v = getValue(field);
        if (v == MISSING_FIELD)
            return defaultValue;
        
        return Integer.parseInt(v);
    }
    
    /**
     * Returns the value of the given field interpreted as Double.
     * <p> Internally the String get converted via 
     * {@link Double#parseDouble(java.lang.String)}.</p>
     * @param field
     * @return 
     */    
    public double getDouble(Integer field) {
        if (getValue(field) == MISSING_FIELD)
            return Double.NaN;
        try {
            return Double.parseDouble(m_currentValues[field]);
        } catch (NumberFormatException nfe) {
            return Double.NaN;
        }
    }

    /**
     * Returns the value of the given field interpreted as Double.
     * <p> Internally the String get converted via 
     * {@link Double#parseDouble(java.lang.String)}.
     * @param field
     * @param defaultValue  if the field is missing this is returned
     * @return 
     */    
    public double getDouble(Integer field, Object defaultValue) {
        if (getValue(field) == MISSING_FIELD)
            return (Double)defaultValue;
        try {
            return Double.parseDouble(m_currentValues[field]);
        } catch (NumberFormatException nfe) {
            return (Double)defaultValue;
        }
    }

   /**
     * Returns the value of the given field interpreted as Boolean.
     * @param field
     * @return null if the field is missing; false if the String value can 
     * be matched to {@link #ISFALSE}; true otherwise
     */    
    public Boolean getBool(Integer field) {
        String v = getValue(field);
        if (v == MISSING_FIELD)
            return null;
        
        return ! ISFALSE.matcher(v).matches();
    }

   /**
     * Returns the value of the given field interpreted as Boolean.
     * @param field
     * @param defaultValue  if the defaultValue is true than anything that can 
     * not be matched to {@link #ISFALSE} will be interpreted as true. Otherwise 
     * only values that can be matched to {@link #ISTRUE} will be result in a 
     * true value;
     * @return null if the field is missing; otherwise according to defaultValue
     * be matched to {@link #ISFALSE}; true otherwise
     */        
    public Boolean getBool(Integer field, boolean defaultValue) {
        String v = getValue(field);
        if (v == MISSING_FIELD)
            return defaultValue;
        
        if (defaultValue)
            return !ISFALSE.matcher(v).matches();
        else
            return ISTRUE.matcher(v).matches();
    }


    /**
     * Returns the value of the given field interpreted as Double.
     * <p> Internally the String get converted via 
     * {@link Double#parseDouble(java.lang.String)}.
     * @param fieldName
     * @return 
     */    
    public double getDouble(String fieldName) {
        int field = getHeaderToColumn().get(fieldName);
        return getDouble(field);
    }

    /**
     * Returns the value of the given field interpreted as Double.
     * <p> Internally the String get converted via 
     * {@link Double#parseDouble(java.lang.String)}.
     * @param fieldName
     * @return 
     */    
    public double getInteger(String fieldName) {
        int field = getHeaderToColumn().get(fieldName);
        return getInteger(field);
    }
    
   /**
     * Returns the value of the given field interpreted as Boolean.
     * @param fieldName
     * @return null if the field is missing; false if the String value can 
     * be matched to {@link #ISFALSE}; true otherwise
     */    
    public Boolean getBool(String fieldName) {
        int field = getHeaderToColumn().get(fieldName);
        return getBool(field);
    }

   /**
     * Returns the value of the given field interpreted as Boolean.
     * @param fieldName
     * @param defaultValue  if the defaultValue is true than anything that can 
     * not be matched to {@link #ISFALSE} will be interpreted as true. Otherwise 
     * only values that can be matched to {@link #ISTRUE} will be result in a 
     * true value;
     * @return null if the field is missing; otherwise according to defaultValue
     * be matched to {@link #ISFALSE}; true otherwise
     */        
    public Boolean getBool(String fieldName, boolean defaultValue) {
        int field = getHeaderToColumn().get(fieldName);
        return getBool(field,defaultValue);
    }
    
    
//    public String[] getFoundValues() {
//        String fv[];
//        fv = java.util.Arrays.<String>copyOf(m_currentValues, m_foundValues);
//        return fv;
//    }
    /**
     * returns an String-array of all values of the current row. 
     * @return 
     */
    public String[] getValues() {
        return m_currentValues;
    }
    
    /**
     * Returns how many values where in the current row.
     * <p>actually it returns number of delimiters + 1.</p>
     * @return 
     */
    public int getFoundValuesCount() {
        return m_foundValues;
    }
    
    /**
     * Returns the highest number of values found in any row of the file so far.
     * <p>This includes the header row</p>
     * @return 
     */
    public int getMaxColumns() {
        return m_maxColumns;
    }

    
    /**
     * Sets the highest number of assumed fields in any row
     * @param columns 
     */
    protected void setMaxColumns(int columns) {
        m_maxColumns=columns;
    }    
    
    /**
     * transfers all information from this CSVParser to the given CSVParser
     * @param csv 
     */
    protected void transfer(CsvParser csv) {
        csv.m_cellValue = m_cellValue;
        csv.m_countHeader = m_countHeader;
        csv.m_currentLine = m_currentLine;
        csv.m_currentValues = m_currentValues.clone();
        csv.m_delimiter = m_delimiter;
        csv.m_doublequote  = m_doublequote;
        csv.m_foundValues = m_foundValues;
        csv.m_header  = m_header;
        csv.m_headerAlternatives  =   (ArrayList<HashSet<String>>) getHeaderAlternatives().clone();
        csv.m_headerToColumn = (HashMap<String,Integer>)getHeaderToColumn().clone();
        csv.m_input = m_input;
        csv.m_inputFile = getInputFile();
        csv.m_maxColumns = m_maxColumns;
        csv.m_quote = m_quote;
        
        
        
    }
    
    /**
     * Returns, whether a header with column-names is defined.
     * @return 
     */
    public boolean hasHeader() {
        return m_header != null;
    }

    /**
     * @return the current row as it was read from the file
     */
    public String getCurrentLine() {
        return m_currentLine;
    }

    /**
     * Defines the givven array as the values for the current row
     * @param currentValues the m_currentValues to set
     */
    protected void setCurrentValues(String[] currentValues) {
        this.m_currentValues = currentValues;
    }
    
    /**
     * Returns the value ready for writing out into a CSV file with the current 
     * configuration.
     * <p>If the value does not contain the quote or delimiter character the 
     * value itself is returned. Otherwise a new String starting and ending 
     * with the quote character is returned. Any occurrence of the quote 
     * character with in the value will be doubled - e.g. <br/> a'b,cd <br/> 
     * would be turned into <br/>'a''b,cd' <br/>. If quote is ' and delimiter is 
     * ,. 
     * @param value
     * @return 
     */
    public String quoteValue(String value) {
        if (value == null)
            return MISSING_FIELD;
        if (value.contains(getDelimiter()) || value.startsWith(getQuote())) {
            return getQuote() + value.replace(getQuote(), getQuote() + getQuote()) + getQuote();
        }
        return value;
    }

    /**
     * Returns the current input file.
     * @return the m_inputFile
     */
    public File getInputFile() {
        return m_inputFile;
    }
    
    /**
     * closes the underlying {@link BufferedReader}.
     * @throws IOException 
     */
    public void close() throws IOException {
        m_input.close();
    }

    /**
     * list of list, containing alternatives for column-names
     * @return the m_headerAlternatives
     */
    public ArrayList<HashSet<String>> getHeaderAlternatives() {
        return m_headerAlternatives;
    }

    /**
     * maps the column-name to the column-index
     * @return the m_headerToColumn
     */
    public HashMap<String,Integer> getHeaderToColumn() {
        return m_headerToColumn;
    }
}
