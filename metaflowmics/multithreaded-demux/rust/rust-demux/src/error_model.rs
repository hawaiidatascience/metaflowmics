// Steps:
// 1) Load barcodes
// 2) hamming_distance(s1, s2)
// 3) for read in fastq:
//      a) get barcodes with dist < min_err
//      b) get all transitions
//      c) add to transition matrix

use itertools::izip;
use seq_io::fastq::{Reader,Record};

extern crate ndarray;
use ndarray::{Array,Array1,Array3,s};

extern crate pav_regression;
pub use pav_regression::pav::{Point,IsotonicRegression};


fn hamming(x: &[u8], y: &[u8]) -> usize {
	x.iter().zip(y.iter())
		.map(|(xi, yi)| (xi != yi) as usize)
		.sum()
}

fn dna_to_int(dna: &u8) -> usize {
	match dna {
		65 => 0,
		67 => 1,
		71 => 2,
		84 => 3,
		78 => 4,
		x => panic!("Unknown character for DNA: {}", x)
	}
}

pub fn count_transitions(
	fastq: &str, barcodes: &Vec<Vec<u8>>, max_dist: usize, n_bases: usize)
	-> Array3<u32>
{
	let mut counts = Array::zeros((5, 4, 41));
	let mut reader = Reader::from_path(fastq).unwrap();
	let mut bp_count = 0;
	
	while bp_count < n_bases {
		match reader.next() {
			None => break,
			Some(read) => {
				let read = read.expect("Error reading record");
				let seq = read.seq();				
				let quals = read.qual();

				for bc in barcodes.iter().filter(|bc| hamming(bc, seq) <= max_dist) {
					for (x, y, q) in izip!(bc, seq, quals) {
						let q = 40.min(*q-33) as usize;
						counts[[dna_to_int(y), dna_to_int(x), q]] += 1;
						bp_count += 1;
					}
				}
			}
		}
	}

	println!("Used {} bp for error model", bp_count);

	// println!("{:?}", counts.slice(s![0, 0, ..]));	
	// println!("{:?}", counts.slice(s![0, 1, ..]));	
	// println!("{:?}", counts.slice(s![0, 2, ..]));
	// println!("{:?}", counts.slice(s![0, 3, ..]));		
	counts
}

pub fn normalize_counts(counts: &mut Array3<f64>) {
	
	// let mut freqs = Array::zeros(counts.raw_dim());

	for i2 in 0..4 {
		for q in 0..41 {
			let new_val: Array1<f64> = {
				let v = counts.slice(s![.., i2, q]);
				let sum = v.sum().max(1.) as f64;
				v.iter().map(|x| *x as f64 / sum).collect()
			};
			// let v = counts.slice(s![.., i2, q]);
			// let sum = v.sum().max(1) as f64;
			// let new_val: Array1<f64> = v.iter().map(|x| *x as f64 / sum).collect();
			counts.slice_mut(s![.., i2, q]).assign(&new_val);
		}
	}
}

pub fn interpolate_counts(counts: &Array3<f64>) -> Array3<f64>{
	let mut interp = Array::zeros(counts.raw_dim());

	for i1 in 0..5 {
		for i2 in 0..4 {
			let x = 0..41;
			let y = counts.slice(s![i1, i2, ..]);

			let points = x.clone().zip(y.iter())
				.filter(|(_, yi)| **yi > 0.)
				.map(|(xi, yi)| Point::new(xi as f64, (*yi as f64).log(10.)))
				.collect::<Vec<Point>>();

			if points.len() < 3 {
				// if i1 < 4 {
				// 	panic!("Too few transitions for NUCL #{} -> NUCL#{}", i1, i2);
				// }
				// continue
				continue
			}

			let model =
				if i1 == i2 {
					IsotonicRegression::new_ascending(&points)
				} else {
					IsotonicRegression::new_descending(&points)
				};
			let y_hat = x
				.map(|xi| model.interpolate(xi as f64).min(0.))
				.map(|xi| (10 as f64).powf(xi))
				.collect::<Array1<f64>>();

			interp.slice_mut(s![i1, i2, ..]).assign(&y_hat);
		}
	}
	interp
}

#[cfg(test)]
mod tests {
	use super::*;

    #[test]
    fn test_hamming() {
		let s1 = "AAAAA".as_bytes().to_vec();
		let s2 = "AATTA".as_bytes().to_vec();
        assert_eq!(hamming(&s1, &s2), 2);
    }
}
