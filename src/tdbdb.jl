function search_db(str::AbstractString)
    strs = split(str)

    N = length(strs)
    elems = String[]
    for i in 1:N
        str = strs[i]
        try
            PeriodicTable.elements[Symbol(str)]
        catch
            continue
        end
        push!(elems, str)
    end
    kwrds = filter(str -> !in(str, elems), strs)

    elems = join(elems, ',')
    res = HTTP.get("http://avdwgroup.engin.brown.edu/getdbid.php?element=" * elems)
    dicts = JSON.parse(String(res.body))

    for kwrd in kwrds
        dicts = filter(dict -> occursin(kwrd, string(dict)), dicts)
    end

    N = length(dicts)
    for i in 1:N
        print("[$i] ")
        dict = dicts[i]

        auth = get(dict, "authoryear", "")
        print(auth * " ")

        elems = get(dict, "element", "")
        elems = join(elems, ',')
        print(elems * " ")

        print('\n')
    end

    if length(dicts) > 1
        print("Select a database to read: ")
        str = readline()
        i = parse(Int, str)
    else
        i = 1
    end

    url = get(dicts[i], "tdburl", "")
    zip = HTTP.download(url)

    archive = ZipFile.Reader(zip)
    @assert length(archive.files) == 1
    file = archive.files[1]
    db = read_tdb(file)
    close(archive)
    Base.Filesystem.rm(zip)

    return db
end