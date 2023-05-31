$1 ~ search {
    for (i = 3; i <= NF; ++i) {
        if ($i == prop) {
            val = $(i + 1)
            gsub(/"/, "", val)
            print val
            exit
        }
    }
}
