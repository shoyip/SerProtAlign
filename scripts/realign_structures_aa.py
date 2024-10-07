if __name__ == '__main__':
    with open('data/3di_aligned.faa', 'r') as al:
        with open('data/3di_aa.faa', 'r') as aa:
            al_lines = al.readlines()
            aa_lines = aa.readlines()

            for al_line, aa_line in zip(al_lines, aa_lines):
                aa_al = ""
                if al_line[0] == ">":
                    print(al_line)
                else:
                    aa_counter = 0
                    for char in al_line:
                        if char == '-':
                            aa_al += '-'
                        else:
                            aa_al += aa_line[aa_counter]
                            aa_counter += 1
                    print(aa_al)
