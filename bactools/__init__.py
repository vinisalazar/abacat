name = "bactools"

from bactools.genome import Genome, load_from_fasta
from bactools.bactools_helper import (
    get_records,
    is_fasta,
    is_fasta_wrapper,
    timer_wrapper,
)
from bactools.prodigal import Prodigal, run
from bactools.prokka import prokka
from bactools.config import CONFIG