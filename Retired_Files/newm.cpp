

#include	"File.h"






int main(int argc, char **argv) {

	int low =  atoi(argv[1]);
	int high = atoi(argv[2]);
	const char *ldir = "temp";

	for(int li = low; li <= high; li++) {

		// and the make file that makes the lists
		//of correspondence points

		char fname[256];
		sprintf( fname, "%s/%d/make.pts", ldir, li );

		FILE	*fp = FileOpenOrDie( fname, "w" );

		fprintf(fp, "all: pts.down pts.same pts.up pts.same.nf\n");
		fprintf(fp, "\n");
		fprintf(fp, "pts.same: */%d*\n", li);
		fprintf(fp, "\talign -fpts.same make.same\n");
		fprintf(fp, "\n");
		fprintf(fp, "pts.same.nf: */%d*\n", li);
		fprintf(fp, "\talign -fpts.same.nf make.same -nf\n");
		fprintf(fp, "\n");
		fprintf(fp, "pts.up: */%d*\n", li+1);
		fprintf(fp, "\talign -fpts.up make.up\n");
		fprintf(fp, "\n");
		fprintf(fp, "pts.down: */%d*\n",li-1);
		fprintf(fp, "\talign -fpts.down make.down\n");
		fclose(fp);
	}

	return 0;
}
