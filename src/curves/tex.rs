pub trait ToTeX<W>
    where W: std::io::Write
{
    fn to_tex(&self, w: &mut W) -> std::io::Result<()>;
}
