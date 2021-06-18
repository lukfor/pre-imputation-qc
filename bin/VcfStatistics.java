import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;

import genepi.io.text.LineWriter;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.variant.vcf.VCFRecordCodec;
import picocli.CommandLine.Help.Visibility;
import picocli.CommandLine;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

//usr/bin/env jbang "$0" "$@" ; exit $?
//REPOS jcenter,jfrog-genepi-maven=https://genepi.jfrog.io/artifactory/maven
//DEPS info.picocli:picocli:4.5.0
//DEPS genepi:genepi-io:1.0.12
//DEPS com.github.samtools:htsjdk:2.21.3

public class VcfStatistics implements Callable<Integer> {

	@Option(names = "--input", description = "input", required = true)
	private String input;

	@Option(names = "--name", description = "name", required = true)
	private String name;

	@Option(names = "--output", description = "output", required = true)
	private String output;

	public static void main(String... args) {
		int exitCode = new CommandLine(new VcfStatistics()).execute(args);
		System.exit(exitCode);
	}

	@Override
	public Integer call() throws Exception {

		assert (output != null);

		int totalSnps = 0;

		VCFFileReader reader = new VCFFileReader(new File(input), false);
		int samples = reader.getFileHeader().getNGenotypeSamples();
		for (VariantContext snp : reader) {
			totalSnps++;
		}

		reader.close();

		// append to existing file
		boolean newFile = false;
		if (!new File(output).exists()) {
			newFile = true;
		}
		BufferedWriter writer = new BufferedWriter(new FileWriter(output, true));
		if (newFile) {
			writer.write("name samples snps");
		}
		writer.write("\n" + name + " " + samples + " " + totalSnps);
		writer.close();

		return 0;
	}

}
