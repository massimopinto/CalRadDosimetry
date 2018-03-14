replace_source='core'
replace_target='strip'
for filename in ./*.dat; do
	new_filename=${filename//$replace_source/$replace_target}
	cut -f 1-2 "$filename" > "$new_filename"
	mv "$new_filename" ./Stripped_data/
done


