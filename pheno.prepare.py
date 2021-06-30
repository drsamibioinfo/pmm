import os
import csv

COLUMNS_DEF = {
    "sample": 0,
    "age": 5,
    "gender": 4,
    "severity": 7,
    "deceased": 8,
    "support": 26,
    "hospitalization": 29
}


def should_exclude(record):
    cols = ["age", "gender", "severity"]
    tt = [len(record[col]) < 1 for col in cols]
    return all(tt)


def get_age(age):
    try:
        int(age)
        return age
    except Exception:
        if "m" in age:
            return 0
        elif "na" in age or "n/a" in age:
            return "n/a"



def get_records(lines):
    for line in lines:
        record = {k: line[v] for k, v in COLUMNS_DEF.items()}
        if not should_exclude(record):
            yield record


def main():
    try:
        input_file = "/media/genomicslab/HD#2/pmm/illumina.tsv"
        out_file = "/media/genomicslab/HD#2/pmm/pmm.pheno.tsv"
        out_headers = ["sample","age","gender","grading","severity"]
        writer = open(out_file,"w")
        writer.write("\t".join(out_headers)+"\n")
        with open(input_file, "r") as input:
            contents = input.readlines()
        contents = [line.replace("\n", "") for line in contents]
        contents = contents[1:]
        lines = [line.split("\t") for line in contents]
        records = list(get_records(lines))
        for record in records:
            if record['sample'] == '772':
                print("Yeah")
            support = record["support"]
            hospitalization = record["hospitalization"]
            hospitalization = hospitalization.lower()
            deceased = record["deceased"].lower()
            support = support.replace(" ","").lower()
            severity = record["severity"].replace(" ","")
            ###################################################### Grading ##############################################
            grading=0
            if hospitalization == "hospitalized":
                grading = 1
            if "non-invasive" in support or "noninvasive" in support or "noneinvasive" in support:
                grading = 2
            elif "invasivemechanical" in support  or "invasiveventilation" in support:
                grading = 4
            if "bibap" in support or "hfnc" in support:
                grading = 3
            if deceased == "1":
                grading = 5
            if len(support) < 1:
                if severity == 'A':
                    grading = 4
                elif severity == 'B':
                    grading = 1
                else:
                    grading = 0
            ################################################### Gender ##################################################
            gt = record["gender"].lower()
            if gt == 'female':
                gender = 1
            else:
                gender = 0
            #################################################### Age ####################################################
            at = record["age"].replace(" ","").lower()
            age = get_age(at)

            onerecord = "\t".join([record["sample"],str(age),str(gender),str(grading),severity])+"\n"
            writer.write(onerecord)
        print("All Finished.")

    except Exception as e:
        print(str(e))


if __name__ == '__main__':
    main()
