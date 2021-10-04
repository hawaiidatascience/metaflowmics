use std::error::Error;

use csv::ReaderBuilder;
use ndarray_npy::write_npy;

use clap::{Arg, App, ArgMatches};

mod error_model;
pub use error_model::{count_transitions,normalize_counts,interpolate_counts};
pub use rust_demux::Barcode;


fn main() {
	let args = parse_args();

	let barcodes = parse_barcodes(
		args.value_of("barcodes").expect("Could not find barcode file")
	).unwrap();
	let fwd = barcodes.into_iter().map(|bc| bc.rev).collect();

	let counts = count_transitions(
		args.value_of("fwd").expect("Could not find fwd fastq file"),
		&fwd, 4, 10000000
	);
	let mut freqs = counts.mapv(|elem| elem as f64);
	normalize_counts(&mut freqs);

	let mut interp = interpolate_counts(&freqs);
	write_npy("freqs_raw.npy", &freqs)
		.expect("Could not write numpy array");

	normalize_counts(&mut interp);
	write_npy("freqs_interp.npy", &interp)
		.expect("Could not write numpy array");

	// println!("{:?}", freqs);
}

fn parse_args() -> ArgMatches<'static> {
    let matches = App::new("rust_demux")
        .version("0.1")
        .about("Read demultiplexing")
        .arg(Arg::with_name("fwd")
             .short("f")
			 .long("fwd")
			 .value_name("FWD")
             .help("Forward fastq")
			 .default_value("data/large/test_I2.fastq")
             .required(true))
        .arg(Arg::with_name("barcodes")
             .short("b")
			 .long("barcodes")
			 .value_name("BARCODE")
             .help("Barcode csv file")
			 .default_value("data/large/sample_barcodes.csv")
             .required(true))
        .get_matches();
    matches
}

fn parse_barcodes(fname: &str) -> Result<Vec<Barcode>, Box<dyn Error>> {
	let mut rdr = ReaderBuilder::new()
        .has_headers(false)
		.delimiter(b',')
		.from_path(fname)?;
	
	let barcodes = rdr
		.records()
		.map(|v| v.unwrap())
		.map(|v| Barcode{
			sample: v[0].to_string(),
			fwd: v[1].as_bytes().to_vec(),
			rev: v[2].as_bytes().to_vec()
		})
		.collect();

	Ok(barcodes)
}
