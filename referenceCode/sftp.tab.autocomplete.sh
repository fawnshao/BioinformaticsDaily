_sftp()
{
    local configfile cur userhost path prefix

    COMPREPLY=()
    cur=`_get_cword ":"`

    _expand || return 0

    if [[ "$cur" == *:* ]]; then
        local IFS=$'\t\n'
        # remove backslash escape from :
        cur=${cur/\\:/:}
        userhost=${cur%%?(\\):*}
        path=${cur#*:}
        # unescape spaces
        path=${path//\\\\\\\\ / }
        if [ -z "$path" ]; then
            # default to home dir of specified user on remote host
            path=$(ssh -o 'Batchmode yes' $userhost pwd 2>/dev/null)
        fi
        # escape spaces; remove executables, aliases, pipes and sockets;
        # add space at end of file names
        COMPREPLY=( $( ssh -o 'Batchmode yes' $userhost \
            command ls -aF1d "$path*" 2>/dev/null | \
            sed -e "s/[][(){}<>\",:;^&\!$=?\`|\\ ']/\\\\\\\\\\\\&/g" \
            -e 's/[*@|=]$//g' -e 's/[^\/]$/& /g' ) )
        return 0
    fi

    if [[ "$cur" = -F* ]]; then
        cur=${cur#-F}
        prefix=-F
    else
        # Search COMP_WORDS for '-F configfile' or '-Fconfigfile' argument
        set -- "${COMP_WORDS[@]}"
        while [ $# -gt 0 ]; do
            if [ "${1:0:2}" = -F ]; then
                if [ ${#1} -gt 2 ]; then
                    configfile="$(dequote "${1:2}")"
                else
                    shift
                    [ "$1" ] && configfile="$(dequote "$1")"
                fi
                break
            fi
            shift
        done

        [[ "$cur" == */* ]] || _known_hosts_real -c -a -F "$configfile" "$cur"
    fi
    # This approach is used instead of _filedir to get a space appended
    # after local file/dir completions, and $nospace retained for others.
    local IFS=$'\t\n'
    COMPREPLY=( "${COMPREPLY[@]}" $( command ls -aF1d $cur* 2>/dev/null | sed \
        -e "s/[][(){}<>\",:;^&\!$=?\`|\\ ']/\\\\&/g" \
        -e 's/[*@|=]$//g' -e 's/[^\/]$/& /g' -e "s/^/$prefix/") )

    return 0
}
complete -o nospace -F _sftp sftp

_complete_ssh_hosts ()
{
        COMPREPLY=()
        cur="${COMP_WORDS[COMP_CWORD]}"
        comp_ssh_hosts=`cat ~/.ssh/known_hosts | \
                        cut -f 1 -d ' ' | \
                        sed -e s/,.*//g | \
                        grep -v ^# | \
                        uniq | \
                        grep -v "\[" ;
                cat ~/.ssh/config | \
                        grep "^Host " | \
                        awk '{print $2}'
                `
        COMPREPLY=( $(compgen -W "${comp_ssh_hosts}" -- $cur))
        return 0
}
complete -F _complete_ssh_hosts ssh