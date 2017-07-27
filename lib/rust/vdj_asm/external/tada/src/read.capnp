@0xb2e4774472d28aa2;

struct BarcodedReadPair {
	gemGroup @0: Int8;
	readName @1: Data;

	read1Seq  @2: Data;
	read1Qual @3: Data;

	read2Seq  @4: Data;
	read2Qual @5: Data;

 	barcode      @6: Data;
	barcodeQual  @7: Data;

	correctedBarcode @8: Data;
	correctionScore  @9: Data;

	trim      @10:  Data;
	trimQual  @11: Data;
}
