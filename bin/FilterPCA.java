
//usr/bin/env jbang "$0" "$@" ; exit $?
//REPOS jcenter,jfrog-genepi-maven=https://genepi.jfrog.io/artifactory/maven
//DEPS info.picocli:picocli:4.5.0
//DEPS com.github.lukfor:magic-tables:0.3.1

import java.util.concurrent.Callable;

import lukfor.tables.Table;
import lukfor.tables.columns.IBuildValueFunction;
import lukfor.tables.columns.types.DoubleColumn;
import lukfor.tables.columns.types.IntegerColumn;
import lukfor.tables.columns.types.StringColumn;
import lukfor.tables.io.TableBuilder;
import lukfor.tables.io.TableWriter;
import lukfor.tables.io.options.CsvTableOptions;
import lukfor.tables.rows.IRowProcessor;
import lukfor.tables.rows.Row;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Help.Visibility;
import picocli.CommandLine.Option;

@Command(name = "MergeTables", description = "Merge multiple csv files into one file")

/**
 * PCA Outlier Removal: Select the number of principal components that should be
 * involved in this process (default 3) and how many standard deviations in all
 * of these components a patient or subject should be within (default 6).
 * Outliners are removed and written to a file. A html report with plots is
 * created.
 *
 * @author lukas
 *
 */
public class FilterPCA implements Callable<Integer> {

	private static final String INDIV_ID = "indivID";

	private static final String SUPERPOP_ID = "superpopID";

	@Option(names = "--reference-pc", description = "input reference_pc file", required = true)
	String referencePc;

	@Option(names = "--reference-samples", description = "input reference samples file", required = true)
	String referenceSamples;

	@Option(names = "--study-pc", description = "input study_pc file", required = true)
	String studyPc;

	@Option(names = "--population", description = "population", required = true)
	String population;

	@Option(names = "--output", description = "output prefix", required = true)
	String output;

	@Option(names = "--max-sd", description = "sd filter", required = false, showDefaultValue = Visibility.ALWAYS)
	int maxSd = 6;

	@Option(names = "--max-pc", description = "max pc used for filtering", required = false, showDefaultValue = Visibility.ALWAYS)
	int maxPC = 3;

	@Option(names = "--column-pc", description = "pc group", required = false, showDefaultValue = Visibility.ALWAYS)
	String columnPc = "PC";

	public static void main(String... args) {
		int exitCode = new CommandLine(new FilterPCA()).execute(args);
		System.exit(exitCode);
	}

	@Override
	public Integer call() throws Exception {

		Table.disableLog();

		Table referenceSampels = TableBuilder.fromCsvFile(new CsvTableOptions(referenceSamples).withSeparator('\t'));

		final Table reference = TableBuilder.fromCsvFile(new CsvTableOptions(referencePc).withSeparator('\t'));
		reference.merge(referenceSampels, INDIV_ID);
		// make unique list of supeor populations
		referenceSampels.getColumns().select(SUPERPOP_ID);
		referenceSampels.getRows().dropDuplicates();

		final Table study = TableBuilder.fromCsvFile(new CsvTableOptions(studyPc).withSeparator('\t'));
		study.getColumns().append(new IntegerColumn("outlier"), new IBuildValueFunction() {
			@Override
			public Object buildValue(Row row) {
				return 0;
			}
		});

		final Table clusters = new Table("PC");
		clusters.getColumns().append(new StringColumn(SUPERPOP_ID));
		for (int i = 1; i <= maxPC; i++) {
			clusters.getColumns().append(new DoubleColumn(columnPc + i + "_mean"));
			clusters.getColumns().append(new DoubleColumn(columnPc + i + "_sd"));
			clusters.getColumns().append(new DoubleColumn(columnPc + i + "_min"));
			clusters.getColumns().append(new DoubleColumn(columnPc + i + "_max"));
		}

		referenceSampels.forEachRow(new IRowProcessor() {
			@Override
			public void process(Row row) {

				String population = row.getString(SUPERPOP_ID);

				System.out.println("Process population " + population + "...");

				Row newRow = clusters.getRows().append();
				newRow.set(0, population);

				for (int i = 0; i < maxPC; i++) {

					System.out.println("  Process " + columnPc + (i + 1) + "...");

					Table referenceCopy = reference.clone();
					referenceCopy.getRows().selectByRegEx(SUPERPOP_ID, population);

					final String column = columnPc + (i + 1);
					double mean = (double) referenceCopy.getColumn(column).getMean();
					double sd = (double) referenceCopy.getColumn(column).getSd();
					final double min = mean - (maxSd * sd);
					final double max = mean + (maxSd * sd);
					newRow.set(i * 4 + 1, mean);
					newRow.set(i * 4 + 2, sd);
					newRow.set(i * 4 + 3, min);
					newRow.set(i * 4 + 4, max);

					if (population.equals(FilterPCA.this.population)) {

						// detect outliners
						study.getColumns().append(new IntegerColumn(column + "_outlier"), new IBuildValueFunction() {
							@Override
							public Object buildValue(Row row) {
								double value = row.getDouble(column);
								boolean outlier = ((value < min) || (value > max));
								if (outlier) {
									row.set("outlier", 1);
								}
								return outlier ? 1 : 0;
							}
						});
					}

				}

			}
		});

		TableWriter.writeToCsv(clusters, output + ".clusters.txt", '\t');

		//keep outliers
		study.getRows().dropByRegEx("outlier", "0");
		TableWriter.writeToCsv(study, output + ".outliers.txt", '\t');


		return 0;
	}

	public void setReferencePc(String referencePc) {
		this.referencePc = referencePc;
	}

	public void setStudyPc(String studyPc) {
		this.studyPc = studyPc;
	}

	public void setOutput(String output) {
		this.output = output;
	}

	public void setMaxPC(int maxPC) {
		this.maxPC = maxPC;
	}

	public void setMaxSd(int maxSd) {
		this.maxSd = maxSd;
	}

	public void setReferenceSamples(String referenceSamples) {
		this.referenceSamples = referenceSamples;
	}

	public void setColumnPc(String columnPc) {
		this.columnPc = columnPc;
	}

	public void setPopulation(String population) {
		this.population = population;
	}

}
