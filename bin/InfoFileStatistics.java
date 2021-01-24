
//usr/bin/env jbang "$0" "$@" ; exit $?
//REPOS jcenter,bintry-lukfor-maven=https://dl.bintray.com/lukfor/maven
//DEPS info.picocli:picocli:4.5.0
//DEPS com.github.lukfor:magic-tables:0.3.1

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;

import genepi.io.FileUtil;
import genepi.io.table.reader.CsvTableReader;
import lukfor.tables.Table;
import lukfor.tables.columns.IBuildValueFunction;
import lukfor.tables.columns.types.DoubleColumn;
import lukfor.tables.columns.types.StringColumn;
import lukfor.tables.io.TableBuilder;
import lukfor.tables.io.TableWriter;
import lukfor.tables.rows.Row;
import lukfor.tables.rows.TableIndex;
import lukfor.tables.rows.filters.IRowFilter;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Help.Visibility;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

@Command(name = "InfoFileStatistics", description = "Bins Rsq by MAF for imputed variants")
public class InfoFileStatistics implements Callable<Integer> {

	public static final double MAX_MAF = 1;

	public static final double MAF_BIN_SIZE = 0.01;

	@Parameters(description = "info files")
	List<String> infoFilenames;

	@Option(names = "--output", description = "the output html filename", required = true)
	String output;

	@Option(names = "--max-maf", description = "filter by MAF", required = false, showDefaultValue = Visibility.ALWAYS)
	double maxMaf = MAX_MAF;

	@Option(names = "--maf-bin-size", description = "MAF bin size", required = false, showDefaultValue = Visibility.ALWAYS)
	double mafBinSize = MAF_BIN_SIZE;

	public static void main(String... args) {
		int exitCode = new CommandLine(new InfoFileStatistics()).execute(args);
		System.exit(exitCode);
	}

	@Override
	public Integer call() throws Exception {

		Table.disableLog();

		Table all = null;


		for (int i = 0; i < infoFilenames.size(); i++) {

			String info = infoFilenames.get(i);
			System.out.println("Read file '" + info + "' (" + (i + 1) + "/" + infoFilenames.size() + ")");

			FileInputStream inputStream = new FileInputStream(info);
			InputStream in = FileUtil.decompressStream(inputStream);
			CsvTableReader reader = new CsvTableReader(new DataInputStream(in), '\t');

			Table table = TableBuilder.fromTableReader("info", reader, false, new String[] { "MAF", "Rsq" });
			table.getRows().dropByRegEx("Rsq", "\\-");
			table.detectTypes();

			// remove all with MAF > maxMaf
			table.getRows().drop(new IRowFilter() {
				@Override
				public boolean accepts(Row row) {
					return row.getDouble("MAF") > maxMaf;
				}
			});

			// bin by MAF and calc mean
			Table binnedMean = table.binBy("MAF", mafBinSize).mean("Rsq");
			binnedMean.getRows().sortAscBy("MAF");
			binnedMean.getColumns().rename("mean", "mean_Rsq");

			// bin by MAF and count values
			Table binnedCount = table.binBy("MAF", mafBinSize).count();
			binnedCount.getRows().sortAscBy("MAF");
			binnedCount.getColumns().rename("count", "count_Rsq");

			binnedMean.merge(binnedCount, "MAF");
			String name = new File(info).getName();
			binnedMean.getColumns().append(new StringColumn("source"), new IBuildValueFunction() {

				@Override
				public Object buildValue(Row row) {
					return name;
				}
			});

			if (all == null) {
				all = binnedMean;
			} else {
				all.append(binnedMean);
			}

		}

		System.out.println("Write result ot file '" + output + "'....");
		TableWriter.writeToCsv(all, output);
		System.out.println("Done.");
		return 0;

	}

	public void setInfoFilenames(String[] infoFilenames) {
		this.infoFilenames = new Vector<String>();
		for (String info : infoFilenames) {
			this.infoFilenames.add(info);
		}
	}

	public void setInfoFilenames(List<String> infoFilenames) {
		this.infoFilenames = infoFilenames;
	}

	public void setMafBinSize(double mafBinSize) {
		this.mafBinSize = mafBinSize;
	}

	public void setMaxMaf(double maxMaf) {
		this.maxMaf = maxMaf;
	}

	public void setOutput(String output) {
		this.output = output;
	}

}
