for ext in cc h py; do 
    find . -type f -name "*$ext" | xargs sed -i "s?${1}?logger.${2}?g"
done
