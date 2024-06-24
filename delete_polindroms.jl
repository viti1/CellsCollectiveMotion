function delete_polindroms( cline::Array{Int,1} )

        
        prev = [length(cline);1:length(cline)-1];
        next = [2:length(cline);1];
        i=1;
        while i<=length(cline)
            #check if polindrom
            prev_j = prev[i];
            next_j = next[i];
            polindrom_size=0;
            while cline[prev_j]==cline[next_j] && polindrom_size<10
                 prev_j = prev[prev_j];
                 next_j = next[next_j];			
                 polindrom_size +=1
            end

            #delete polindrom
            if polindrom_size>0
				println("polindrom_size= $polindrom_size")
                if i+polindrom_size >length(cline)
					range1 = i-polindrom_size+1:length(cline);
					range2 = 1:i+polindrom_size-length(cline);
                    splice!(cline,range1)
                    splice!(cline,range2)
					break;
                elseif i-polindrom_size+1 < 1
					range1 = length(cline)+i-polindrom_size+1:length(cline);
					range2 = 1:i+polindrom_size;
					splice!(cline,range1)
					splice!(cline,range2)
					i=1
				else
					splice!(cline,i-polindrom_size+1:i+polindrom_size)
					i=i-polindrom_size+1;
				end
                prev = [length(cline);1:length(cline)-1];
                next = [2:length(cline);1];
            else
				i+=1;
			end
        end
		
#		println(cline')
end
