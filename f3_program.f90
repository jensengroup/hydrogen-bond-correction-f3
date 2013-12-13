program thirdgeneration_hbond_correction
    ! """
    ! Third Generation Hydrogen-bond Correction
    ! https://github.com/jensengroup/hydrogen-bond-correction-f3
    ! """

    implicit none

    ! Version
    character(len=*), parameter :: version = '1.0'

    ! Periodic table
    character(len=2) :: element(94)
    data element /'h ','he', &
    & 'li','be','b ','c ','n ','o ','f ','ne', &
    & 'na','mg','al','si','p ','s ','cl','ar', &
    & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
    & 'zn','ga','ge','as','se','br','kr', &
    & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
    & 'cd','in','sn','sb','te','i ','xe', &
    & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
    & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
    & 'au','hg','tl','pb','bi','po','at','rn', &
    & 'fr','ra','ac','th','pa','u ','np','pu'/

    ! Options
    character(len=32) :: arg
    character(len=32) :: paramfile, geofile
    character(len=70) :: line
    logical :: param, debug, lgradient
    logical :: exists, skip = .False.
    integer :: success, ix, ierr

    ! Model
    integer :: natoms
    integer, parameter :: maxnatoms = 2000

    character(len=2) :: labels(maxnatoms)
    double precision :: geo(3,maxnatoms)
    double precision :: gradient(3,maxnatoms)
    integer :: ilabels(maxnatoms)

    double precision :: c_oxygen, c_nitrogen

    integer :: i, j

    param = .False.
    debug = .False.
    lgradient = .False.
    geofile = ''

    if(command_argument_count() == 0) then
        call print_help()
        stop
    endif

    do i = 1, command_argument_count()

        if(skip)then
            skip=.False.
            cycle
        endif

        call get_command_argument(i, arg)

        select case (arg)
        case ('-g', '--gradient')
            lgradient = .True.

        case ('-d', '--debug')
            debug = .True.

        case ('-p', '--param')
            param = .True.

            call get_command_argument(i+1, paramfile)
            inquire(file=paramfile, exist=exists)
            if(.not.exists) then
                print '(a)', ''
                print '(2a)', 'Parameter file missing ', paramfile
                print '(a)', ''
                stop
            else
                skip = .True.
            endif

        case ('-v', '--version')
            print '(a)', ''
            print '(2a)', 'version ', version
            print '(a)', 'https://github.com/jensengroup/hydrogen-bond-correction-f3'
            call print_cite()
            stop

        case ('-h', '--help')
            call print_help()
            stop

        case default
            inquire(file=arg, exist=exists)
            if(.not.exists) then
                print '(a)', ''
                print '(2a)', 'Unrecognized option or file: ', arg
                print '(a)', ''
                stop
            else
                geofile = arg
            endif

        end select
    enddo


    if(geofile=='')then
        call print_error('No structure file')
    endif

    !
    ! Open Coordinates
    ! Code is optimized for reading XYZ format
    !
    open(ix, file=geofile, iostat=ierr)
    if(ierr.ne.0)then
        call print_error('Unable to open input file')
    endif

    read(ix,*) natoms
    if(natoms.gt.maxnatoms) then
        call print_error('Too many atoms')
    endif

    read(ix,*) ! I'm not interested in the name

    do i=1,natoms

        ! Read XYZ format
        ! a, float, float, float
        read(ix, *, iostat=success) labels(i), geo(1,i), geo(2,i), geo(3,i)
        if (success.ne.0) exit

        ! Convert atom-type from char to int
        call lower_case(labels(i))
        do j = 1, 94
            if(element(j)==labels(i)) then
                ilabels(i) = j
                exit
            endif
        enddo

    enddo

    ! Open Parameters
    if(param)then
        open(ix, file=paramfile, iostat=ierr)
        if(ierr.ne.0)then
            call print_error('Unable to open parameter file')
        endif

        read(ix, *, iostat=success) c_nitrogen
        if (success.ne.0) call print_error('Missing Nitrogen parameter')

        read(ix, *, iostat=success) c_oxygen
        if (success.ne.0) call print_error('Missing Oxygen parameter')
    else
        ! Default parameters
        c_nitrogen = -0.11d0
        c_oxygen   = -0.15d0
    endif

    ! Calculate energy

    ! Calculate gradient

    ! Print output
    print '(a)', ''
    print '(a)', ''
    print '(a)', '*************************************'
    print '(a)', '      Hydrogen-Bond Correction'
    print '(a)', '*************************************'
    print '(a)', ''
    print '(a)', 'Energy:'
    print '(a)', '      -456.32    kcal/mol'
    print '(a)', '      -256.32    joule'
    print '(a)', '      -156.32    hatree'
    print '(a)', ''
    call print_cite()

contains

    subroutine print_error(msg)
        implicit none
        character (len=*) , intent(in) :: msg
        print '(a)', ''
        print '(a)', msg
        print '(a)', ''
        stop
    end subroutine print_error

    subroutine print_cite()
        print '(a)', ''
        print '(a)', 'Cite as:'
        print '(a)', 'Martin Korth'
        print '(a)', 'J. Chem. Theory Comput., 2010, 6 (12), pp 3808–3816'
        print '(a)', 'DOI: 10.1021/ct100408b'
        print '(a)', ''
    end subroutine print_cite

    subroutine print_help()
        print '(a)', ''
        print '(a)', 'usage: f3 [OPTIONS] <structure.xyz>'
        print '(a)', ''
        print '(a)', 'Calculate Hydrogen-bond correction for a XYZ structure'
        print '(a)', ''
        print '(a)', 'options:'
        print '(a)', ''
        print '(a)', '  -v, --version     print version information and exit'
        print '(a)', '  -h, --help        print usage information and exit'
        print '(a)', '  -p <file>, --param <file>   use parameter file <par.dat>'
        print '(a)', '  -d, --debug       print debug information'
        print '(a)', '  -g, --gradient    calculate and print gradient'
        print '(a)', ''
    end subroutine print_help

    subroutine lower_case(word)
        character (len=*) , intent(in out) :: word
        integer :: i, ic, nlen
        nlen = len(word)
        do i = 1, nlen
            ic = ichar(word(i:i))
            if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
        end do
    end subroutine lower_case

end program thirdgeneration_hbond_correction
